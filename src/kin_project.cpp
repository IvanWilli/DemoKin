// Projection kernels for DemoKin.
//
// All kinship models in the package are built from the same two operations:
//   (1) y = U x + F z   (project a kin vector one age/time step, add newborns)
//   (2) k[,1] = M %*% pi (initial condition: average a kin matrix over the
//                         parents' age-at-birth distribution)
// The original R code performs (1) with dense `%*%` on matrices that are
// >99% zeros, and R's `%*%` (reference BLAS dgemv) accumulates each output
// element in ascending column order, starting from 0. The kernels below do
// exactly that accumulation over the non-zero entries only, so results are
// bit-identical to the dense products while costing O(nnz) instead of O(n^2).
//
// Floating-point contract: fused multiply-add is disabled below (pragmas,
// equivalent to -ffp-contract=off, which R CMD check would reject as a
// non-portable Makevars flag) so that no FMA changes the rounding of a*x + s.

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <cstring>

using namespace Rcpp;

namespace {

// Compressed sparse row storage with column indices ascending within a row.
struct CSR {
  int nrow = 0;
  int ncol = 0;
  std::vector<int> ptr;
  std::vector<int> idx;
  std::vector<double> val;
  std::vector<int> nzrows;   // rows with at least one entry, ascending
  void index_rows() {
    nzrows.clear();
    for (int r = 0; r < nrow; ++r) if (ptr[r + 1] > ptr[r]) nzrows.push_back(r);
  }
};

// Build CSR from a dense matrix: every non-zero entry is kept, in column order.
// One pass over the matrix collects the non-zeros (column-major, so ascending
// column within each row) and counts them per row; a counting sort then lays
// them out by row, keeping that order.
CSR csr_from_dense(const NumericMatrix& A) {
  CSR S;
  S.nrow = A.nrow();
  S.ncol = A.ncol();
  S.ptr.assign(S.nrow + 1, 0);
  const double* a = A.begin();
  std::vector<int> ri, ci;
  std::vector<double> vv;
  for (int j = 0; j < S.ncol; ++j) {
    const double* aj = a + (size_t)j * S.nrow;
    for (int i = 0; i < S.nrow; ++i) {
      const double v = aj[i];
      if (v != 0.0) { ri.push_back(i); ci.push_back(j); vv.push_back(v); S.ptr[i + 1]++; }
    }
  }
  for (int i = 0; i < S.nrow; ++i) S.ptr[i + 1] += S.ptr[i];
  S.idx.resize(S.ptr[S.nrow]);
  S.val.resize(S.ptr[S.nrow]);
  std::vector<int> fill(S.ptr.begin(), S.ptr.end() - 1);
  for (size_t k = 0; k < ri.size(); ++k) {
    const int p = fill[ri[k]]++;
    S.idx[p] = ci[k];
    S.val[p] = vv[k];
  }
  S.index_rows();
  return S;
}

// Build CSR from triplets (1-based i, j from R). Entries are sorted by
// (row, column); duplicates are not expected and are summed in column order.
CSR csr_from_triplets(const IntegerVector& i, const IntegerVector& j,
                      const NumericVector& x, int nrow, int ncol) {
  CSR S;
  S.nrow = nrow;
  S.ncol = ncol;
  const int nnz = i.size();
  std::vector<int> ord(nnz);
  for (int k = 0; k < nnz; ++k) ord[k] = k;
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (i[a] != i[b]) return i[a] < i[b];
    return j[a] < j[b];
  });
  S.ptr.assign(nrow + 1, 0);
  S.idx.reserve(nnz);
  S.val.reserve(nnz);
  for (int k = 0; k < nnz; ++k) {
    int o = ord[k];
    if (x[o] == 0.0) continue;
    S.ptr[i[o]]++;  // i is 1-based, so this counts into row i[o]-1's slot + 1
    S.idx.push_back(j[o] - 1);
    S.val.push_back(x[o]);
  }
  for (int r = 0; r < nrow; ++r) S.ptr[r + 1] += S.ptr[r];
  S.index_rows();
  return S;
}

// Accept either a dense numeric matrix or a list(i, j, x, nrow, ncol).
CSR csr_from_sexp(SEXP A) {
  if (Rf_isMatrix(A)) return csr_from_dense(NumericMatrix(A));
  List L(A);
  return csr_from_triplets(L["i"], L["j"], L["x"], as<int>(L["nrow"]), as<int>(L["ncol"]));
}

// y = A x, accumulating in ascending column order from 0.0 (reference dgemv).
inline void matvec(const CSR& A, const double* x, double* y) {
  for (int r = 0; r < A.nrow; ++r) {
    double s = 0.0;
    for (int k = A.ptr[r]; k < A.ptr[r + 1]; ++k) s += A.val[k] * x[A.idx[k]];
    y[r] = s;
  }
}

// y = A x + B z  (two dgemv calls added elementwise, as in R)
// Rows where B has no entry contribute (B z)[r] = 0.0, and y[r] + 0.0 is
// y[r] bit for bit: y[r] is a sum started at +0.0, so it is never -0.0, and
// NaN and Inf are unchanged by adding 0. Only B's non-empty rows are visited.
inline void matvec2(const CSR& A, const double* x, const CSR& B, const double* z,
                    double* y, std::vector<double>& tmp) {
  (void)tmp;
  matvec(A, x, y);
  for (int r : B.nzrows) {
    double s = 0.0;
    for (int k = B.ptr[r]; k < B.ptr[r + 1]; ++k) s += B.val[k] * z[B.idx[k]];
    y[r] = y[r] + s;
  }
}

// out[0:nlive) = M[0:nlive, ] %*% w   (dense dgemv, ascending column order)
inline void ic_product(const NumericMatrix& M, int nlive, const double* w, int nw, double* out) {
  const int n = M.nrow();
  const double* m = M.begin();
  for (int r = 0; r < nlive; ++r) {
    double s = 0.0;
    for (int j = 0; j < nw; ++j) s += m[r + (size_t)j * n] * w[j];
    out[r] = s;
  }
}

inline double* col(NumericMatrix& X, int j) { return X.begin() + (size_t)j * X.nrow(); }
inline const double* ccol(const NumericMatrix& X, int j) { return X.begin() + (size_t)j * X.nrow(); }

}  // namespace

//' Focal's trip through the kinship network (time invariant or one time step)
//'
//' Internal kernel. Replicates the loops of `kin_time_invariant()`,
//' `timevarying_kin()` and their two-sex / cause-of-death / multistage variants.
//'
//' @param U,F,Fstar,G dense matrices or `list(i, j, x, nrow, ncol)` triplets.
//' @param phi0 initial Focal vector (length n).
//' @param m0 initial parents vector placed into `m[seq_along(m0), 1]`.
//' @param pi_all distribution used for gm, ggm, oa, ya, coa, cya initial conditions.
//' @param pi_os distribution used for os, nos initial conditions.
//' @param ncol number of Focal ages (columns).
//' @param nlive number of rows that enter the initial-condition products.
//' @param prev previous period's kin list (time varying step) or `NULL`.
//' @return named list of n x ncol matrices: Focal, d, gd, ggd, m, gm, ggm, os,
//'   ys, nos, nys, oa, ya, coa, cya.
//' @noRd
// [[Rcpp::export]]
List kin_project_cpp(SEXP U, SEXP F, SEXP Fstar, SEXP G,
                     NumericVector phi0, NumericVector m0,
                     NumericVector pi_all, NumericVector pi_os,
                     int ncol, int nlive, Nullable<List> prev) {
  const CSR Uc = csr_from_sexp(U), Fc = csr_from_sexp(F), Fsc = csr_from_sexp(Fstar), Gc = csr_from_sexp(G);
  const int n = Uc.nrow;
  const bool tv = prev.isNotNull();
  const int om = ncol - 1;
  std::vector<double> tmp(n);

  auto zeros = [&]() { return NumericMatrix(n, ncol); };
  NumericMatrix phi = zeros(), d = zeros(), gd = zeros(), ggd = zeros(), m = zeros(), gm = zeros(),
                ggm = zeros(), os = zeros(), ys = zeros(), nos = zeros(), nys = zeros(), oa = zeros(),
                ya = zeros(), coa = zeros(), cya = zeros();

  // sources for the right-hand sides: previous period (time varying) or current
  NumericMatrix pd, pgd, pggd, pm, pgm, pggm, pos, pys, pnos, pnys, poa, pya, pcoa, pcya;
  if (tv) {
    List P(prev);
    pd = as<NumericMatrix>(P["d"]);   pgd = as<NumericMatrix>(P["gd"]);   pggd = as<NumericMatrix>(P["ggd"]);
    pm = as<NumericMatrix>(P["m"]);   pgm = as<NumericMatrix>(P["gm"]);   pggm = as<NumericMatrix>(P["ggm"]);
    pos = as<NumericMatrix>(P["os"]); pys = as<NumericMatrix>(P["ys"]);   pnos = as<NumericMatrix>(P["nos"]);
    pnys = as<NumericMatrix>(P["nys"]); poa = as<NumericMatrix>(P["oa"]); pya = as<NumericMatrix>(P["ya"]);
    pcoa = as<NumericMatrix>(P["coa"]); pcya = as<NumericMatrix>(P["cya"]);
  } else {
    pd = d; pgd = gd; pggd = ggd; pm = m; pgm = gm; pggm = ggm; pos = os; pys = ys;
    pnos = nos; pnys = nys; poa = oa; pya = ya; pcoa = coa; pcya = cya;
  }

  // initial conditions
  for (int r = 0; r < n; ++r) phi(r, 0) = phi0[r];
  for (int r = 0; r < m0.size(); ++r) m(r, 0) = m0[r];
  const int npi = pi_all.size(), nos_pi = pi_os.size();

  if (tv) {
    ic_product(pm, nlive, pi_all.begin(), npi, col(gm, 0));
    ic_product(pgm, nlive, pi_all.begin(), npi, col(ggm, 0));
    ic_product(pos, nlive, pi_all.begin(), npi, col(oa, 0));
    ic_product(pys, nlive, pi_all.begin(), npi, col(ya, 0));
    ic_product(pnos, nlive, pi_all.begin(), npi, col(coa, 0));
    ic_product(pnys, nlive, pi_all.begin(), npi, col(cya, 0));
    ic_product(pd, nlive, pi_os.begin(), nos_pi, col(os, 0));
    ic_product(pgd, nlive, pi_os.begin(), nos_pi, col(nos, 0));
  }

  // first pass: Focal, descendants, parents, younger siblings, their children
  for (int i = 0; i < om; ++i) {
    matvec(Gc, ccol(phi, i), col(phi, i + 1));
    matvec2(Uc, ccol(pd, i), Fc, ccol(phi, i), col(d, i + 1), tmp);
    matvec2(Uc, ccol(pgd, i), Fc, ccol(pd, i), col(gd, i + 1), tmp);
    matvec2(Uc, ccol(pggd, i), Fc, ccol(pgd, i), col(ggd, i + 1), tmp);
    matvec(Uc, ccol(pm, i), col(m, i + 1));
    matvec2(Uc, ccol(pys, i), Fsc, ccol(pm, i), col(ys, i + 1), tmp);
    matvec2(Uc, ccol(pnys, i), Fc, ccol(pys, i), col(nys, i + 1), tmp);
  }

  // grandparents, great-grandparents
  if (!tv) ic_product(m, nlive, pi_all.begin(), npi, col(gm, 0));
  for (int i = 0; i < om; ++i) matvec(Uc, ccol(pgm, i), col(gm, i + 1));
  if (!tv) ic_product(gm, nlive, pi_all.begin(), npi, col(ggm, 0));
  for (int i = 0; i < om; ++i) matvec(Uc, ccol(pggm, i), col(ggm, i + 1));

  // older siblings and their children
  if (!tv) {
    ic_product(d, nlive, pi_os.begin(), nos_pi, col(os, 0));
    ic_product(gd, nlive, pi_os.begin(), nos_pi, col(nos, 0));
  }
  for (int i = 0; i < om; ++i) {
    matvec(Uc, ccol(pos, i), col(os, i + 1));
    matvec2(Uc, ccol(pnos, i), Fc, ccol(pos, i), col(nos, i + 1), tmp);
  }

  // aunts and cousins
  if (!tv) {
    ic_product(os, nlive, pi_all.begin(), npi, col(oa, 0));
    ic_product(ys, nlive, pi_all.begin(), npi, col(ya, 0));
    ic_product(nos, nlive, pi_all.begin(), npi, col(coa, 0));
    ic_product(nys, nlive, pi_all.begin(), npi, col(cya, 0));
  }
  for (int i = 0; i < om; ++i) {
    matvec(Uc, ccol(poa, i), col(oa, i + 1));
    matvec2(Uc, ccol(pya, i), Fsc, ccol(pgm, i), col(ya, i + 1), tmp);
    matvec2(Uc, ccol(pcoa, i), Fc, ccol(poa, i), col(coa, i + 1), tmp);
    matvec2(Uc, ccol(pcya, i), Fc, ccol(pya, i), col(cya, i + 1), tmp);
  }

  return List::create(_["Focal"] = phi, _["d"] = d, _["gd"] = gd, _["ggd"] = ggd, _["m"] = m,
                      _["gm"] = gm, _["ggm"] = ggm, _["os"] = os, _["ys"] = ys, _["nos"] = nos,
                      _["nys"] = nys, _["oa"] = oa, _["ya"] = ya, _["coa"] = coa, _["cya"] = cya);
}

//' Two-sex multistage projection (Butterick), time invariant or one time step
//'
//' Internal kernel replicating `all_kin_dy()` / `all_kin_dy_TV()` exactly,
//' including the order in which initial conditions are accumulated.
//'
//' @param Uproj,Fproj,Fprojstar,Gtilde dense matrices or triplet lists.
//' @param IC_Focal Focal's initial vector.
//' @param parents parents' joint age-stage-sex distribution.
//' @param pif,pim mothers' and fathers' marginal age distributions (length na).
//' @param na number of ages.
//' @param prev previous period's kin list or `NULL` (boundary year).
//' @noRd
// [[Rcpp::export]]
List kin_project_ms2sex_cpp(SEXP Uproj, SEXP Fproj, SEXP Fprojstar, SEXP Gtilde,
                            NumericVector IC_Focal, NumericVector parents,
                            NumericVector pif, NumericVector pim, int na, Nullable<List> prev) {
  const CSR Uc = csr_from_sexp(Uproj), Fc = csr_from_sexp(Fproj), Fsc = csr_from_sexp(Fprojstar), Gc = csr_from_sexp(Gtilde);
  const int n2 = Uc.nrow;
  const bool tv = prev.isNotNull();
  std::vector<double> tmp(n2);

  auto zeros = [&]() { return NumericMatrix(n2, na); };
  NumericMatrix Fo = zeros(), ch = zeros(), gch = zeros(), ggch = zeros(), pa = zeros(), gpa = zeros(),
                ggpa = zeros(), osb = zeros(), ysb = zeros(), onn = zeros(), ynn = zeros(), oau = zeros(),
                yau = zeros(), oco = zeros(), yco = zeros();

  NumericMatrix pFo, pch, pgch, pggch, ppa, pgpa, pggpa, posb, pysb, ponn, pynn, poau, pyau, poco, pyco;
  if (tv) {
    List P(prev);
    pFo = as<NumericMatrix>(P["Focal"]); pch = as<NumericMatrix>(P["d"]); pgch = as<NumericMatrix>(P["gd"]);
    pggch = as<NumericMatrix>(P["ggd"]); ppa = as<NumericMatrix>(P["m"]); pgpa = as<NumericMatrix>(P["gm"]);
    pggpa = as<NumericMatrix>(P["ggm"]); posb = as<NumericMatrix>(P["os"]); pysb = as<NumericMatrix>(P["ys"]);
    ponn = as<NumericMatrix>(P["nos"]); pynn = as<NumericMatrix>(P["nys"]); poau = as<NumericMatrix>(P["oa"]);
    pyau = as<NumericMatrix>(P["ya"]); poco = as<NumericMatrix>(P["coa"]); pyco = as<NumericMatrix>(P["cya"]);
  } else {
    pFo = Fo; pch = ch; pgch = gch; pggch = ggch; ppa = pa; pgpa = gpa; pggpa = ggpa; posb = osb;
    pysb = ysb; ponn = onn; pynn = ynn; poau = oau; pyau = yau; poco = oco; pyco = yco;
  }

  for (int r = 0; r < n2; ++r) { Fo(r, 0) = IC_Focal[r]; pa(r, 0) = parents[r]; }

  // y[,1] = y[,1] + c * X[,ic]
  auto accumulate = [&](NumericMatrix& Y, double c, const NumericMatrix& X, int ic) {
    double* y = col(Y, 0);
    const double* x = ccol(X, ic);
    for (int r = 0; r < n2; ++r) y[r] = y[r] + c * x[r];
  };

  // projection of kin with deterministic initial conditions
  for (int i = 0; i < na - 1; ++i) {
    matvec(Gc, ccol(pFo, i), col(Fo, i + 1));
    matvec(Uc, ccol(ppa, i), col(pa, i + 1));
    matvec2(Uc, ccol(pysb, i), Fsc, ccol(ppa, i), col(ysb, i + 1), tmp);
    matvec2(Uc, ccol(pynn, i), Fc, ccol(pysb, i), col(ynn, i + 1), tmp);
    matvec2(Uc, ccol(pch, i), Fc, ccol(pFo, i), col(ch, i + 1), tmp);
    matvec2(Uc, ccol(pgch, i), Fc, ccol(pch, i), col(gch, i + 1), tmp);
    matvec2(Uc, ccol(pggch, i), Fc, ccol(pgch, i), col(ggch, i + 1), tmp);
  }

  // initial conditions derived from the above (same statement order as the R code:
  // in the boundary year great-grandparents read the partially filled grandparents)
  for (int ic = 0; ic < na; ++ic) {
    double c = pif[ic] + pim[ic];
    accumulate(gpa, c, ppa, ic);
    accumulate(ggpa, c, pgpa, ic);
    accumulate(osb, pif[ic], pch, ic);
    accumulate(onn, pif[ic], pgch, ic);
    accumulate(yau, c, pysb, ic);
  }
  for (int i = 0; i < na - 1; ++i) {
    matvec(Uc, ccol(pgpa, i), col(gpa, i + 1));
    matvec(Uc, ccol(pggpa, i), col(ggpa, i + 1));
    matvec(Uc, ccol(posb, i), col(osb, i + 1));
    matvec2(Uc, ccol(ponn, i), Fc, ccol(posb, i), col(onn, i + 1), tmp);
    matvec2(Uc, ccol(pyau, i), Fsc, ccol(pgpa, i), col(yau, i + 1), tmp);
  }

  // older aunts/uncles and cousins: the R loop runs over 1:(na-1) only
  for (int ic = 0; ic < na - 1; ++ic) {
    double c = pif[ic] + pim[ic];
    accumulate(oau, c, posb, ic);
    accumulate(oco, c, ponn, ic);
    accumulate(yco, c, pynn, ic);
  }
  for (int i = 0; i < na - 1; ++i) {
    matvec(Uc, ccol(poau, i), col(oau, i + 1));
    matvec2(Uc, ccol(poco, i), Fc, ccol(poau, i), col(oco, i + 1), tmp);
    matvec2(Uc, ccol(pyco, i), Fc, ccol(pyau, i), col(yco, i + 1), tmp);
  }

  return List::create(_["Focal"] = Fo, _["d"] = ch, _["gd"] = gch, _["ggd"] = ggch, _["m"] = pa,
                      _["gm"] = gpa, _["ggm"] = ggpa, _["os"] = osb, _["ys"] = ysb, _["nos"] = onn,
                      _["nys"] = ynn, _["oa"] = oau, _["ya"] = yau, _["coa"] = oco, _["cya"] = yco);
}

//' Sum over ages within (stage, Focal age) cells using long double accumulation
//'
//' Replicates `dplyr::summarise(num = sum(value))` (R's `sum()` accumulates in
//' long double) over a block of rows laid out age-major, stage fastest.
//'
//' @param X matrix with `na * ns` rows (stage fastest) and `nc` columns.
//' @param ns number of stages.
//' @return ns x nc matrix of sums.
//' @noRd
// [[Rcpp::export]]
NumericMatrix stage_sums_cpp(NumericMatrix X, int ns) {
  const int nr = X.nrow(), nc = X.ncol();
  const int na = nr / ns;
  NumericMatrix out(ns, nc);
  for (int j = 0; j < nc; ++j) {
    const double* x = ccol(X, j);
    for (int s = 0; s < ns; ++s) {
      long double acc = 0.0;
      for (int a = 0; a < na; ++a) acc += x[a * ns + s];
      out(s, j) = (double)acc;
    }
  }
  return out;
}

//' Grouped kin summary with R's accumulation semantics
//'
//' Replicates the `dplyr::summarise()` / `cumsum()` block of `kin()` exactly:
//' sums accumulate in long double over rows in the given order (as R's `sum()`
//' and `cumsum()` do) and the standard deviation uses `R_pow(x, 0.5)` like `^.5`.
//'
//' @param living,dead,age_kin columns of `kin_full`, in their original row order.
//' @param ord ordering of the rows (1-based) that makes each group contiguous.
//' @param start 1-based positions in `ord` where groups start.
//' @param ngroups number of groups.
//' @param age_focal_g Focal's age of each group.
//' @param cg 0-based id of the cumulative group (kin x sex x period/cohort) of each group.
//' @param ncg number of cumulative groups.
//' @noRd
// [[Rcpp::export]]
List kin_summary_cpp(NumericVector living, NumericVector dead, SEXP age_kin,
                     IntegerVector ord, IntegerVector start, int ngroups, NumericVector age_focal_g,
                     IntegerVector cg, int ncg) {
  // rows are read through ord (1-based), group by group: the same values in
  // the same order as the sorted copies the R side used to pass
  const R_xlen_t n = ord.size();
  const bool age_int = TYPEOF(age_kin) == INTSXP;
  const int* ai = age_int ? INTEGER(age_kin) : nullptr;
  const double* ad = age_int ? nullptr : REAL(age_kin);
  std::vector<long double> s_live(ngroups, 0.0L), s_la(ngroups, 0.0L), s_la2(ngroups, 0.0L), s_dead(ngroups, 0.0L);
  for (int k = 0; k < ngroups; ++k) {
    const R_xlen_t a0 = (R_xlen_t)start[k] - 1;
    const R_xlen_t a1 = (k + 1 < ngroups) ? (R_xlen_t)start[k + 1] - 1 : n;
    for (R_xlen_t p = a0; p < a1; ++p) {
      const R_xlen_t i = (R_xlen_t)ord[p] - 1;
      const double a = age_int ? (ai[i] == NA_INTEGER ? NA_REAL : (double)ai[i]) : ad[i];
      s_live[k] += living[i];
      s_la[k] += living[i] * a;
      s_la2[k] += living[i] * (a * a);
      s_dead[k] += dead[i];
    }
  }
  NumericVector count_living(ngroups), mean_age(ngroups), sd_age(ngroups), count_dead(ngroups),
      count_cum_dead(ngroups), mean_age_lost(ngroups);
  std::vector<long double> cum_dead(ncg, 0.0L), cum_dead_age(ncg, 0.0L);
  for (int k = 0; k < ngroups; ++k) {
    const double cl = (double)s_live[k];
    const double la = (double)s_la[k];
    const double la2 = (double)s_la2[k];
    count_living[k] = cl;
    const double ma = la / cl;
    mean_age[k] = ma;
    sd_age[k] = R_pow(la2 / cl - ma * ma, 0.5);
    const double cd = (double)s_dead[k];
    count_dead[k] = cd;
    const int c = cg[k];
    cum_dead[c] += cd;
    cum_dead_age[c] += cd * age_focal_g[k];
    count_cum_dead[k] = (double)cum_dead[c];
    mean_age_lost[k] = (double)cum_dead_age[c] / (double)cum_dead[c];
  }
  return List::create(_["count_living"] = count_living, _["mean_age"] = mean_age, _["sd_age"] = sd_age,
                      _["count_dead"] = count_dead, _["count_cum_dead"] = count_cum_dead,
                      _["mean_age_lost"] = mean_age_lost);
}

//' Grouped sums with R's `sum()` semantics
//'
//' Replicates `sum(x)` evaluated group by group, as `dplyr::summarise()` does:
//' within a group the values are added in row order into a long double that
//' starts at 0 (R's `rsum()`), the total is clamped to +/-Inf beyond
//' DBL_MAX, rounded to double, and added to 0.0 (the `zcum.r += tmp` step of
//' `do_summary()`). NA and NaN propagate as in R, since the arithmetic is the
//' same.
//'
//' @param x values, in the original row order.
//' @param o ordering of the rows (1-based) that makes each group's rows
//'   contiguous and keeps their original order within a group.
//' @param start 1-based positions in `o` where groups start, increasing.
//' @return one sum per group.
//' @noRd
// [[Rcpp::export]]
NumericVector group_sum_cpp(NumericVector x, IntegerVector o, IntegerVector start) {
  const R_xlen_t n = o.size();
  const R_xlen_t ng = start.size();
  NumericVector out(ng);
  const double* px = x.begin();
  const int* po = o.begin();
  for (R_xlen_t g = 0; g < ng; ++g) {
    const R_xlen_t a = (R_xlen_t)start[g] - 1;
    const R_xlen_t b = (g + 1 < ng) ? (R_xlen_t)start[g + 1] - 1 : n;
    long double s = 0.0;
    for (R_xlen_t i = a; i < b; ++i) s += px[po[i] - 1];
    double v;
    if (s > DBL_MAX) v = R_PosInf;
    else if (s < -DBL_MAX) v = R_NegInf;
    else v = (double)s;
    double z = 0.0;
    z += v;
    out[g] = z;
  }
  return out;
}

//' First row of each group of sorted keys
//'
//' `o` is a (1-based) ordering of the rows that makes every key column
//' non-decreasing (e.g. `order(..., method = "radix")`). Row `o[i]` starts a
//' new group when any key differs from row `o[i-1]`; two values are the same
//' when they are equal or both missing (NA or NaN), which is what grouping
//' does. Strings are compared by content.
//'
//' @param keys list of key columns (integer, logical, double or character).
//' @param o ordering of the rows, 1-based.
//' @return 1-based positions in `o` where groups start.
//' @noRd
// [[Rcpp::export]]
IntegerVector group_starts_cpp(List keys, IntegerVector o) {
  const R_xlen_t n = o.size();
  if (n == 0) return IntegerVector(0);
  std::vector<char> brk(n, 0);
  brk[0] = 1;
  const int* po = o.begin();
  for (R_xlen_t c = 0; c < keys.size(); ++c) {
    SEXP x = keys[c];
    switch (TYPEOF(x)) {
    case INTSXP: case LGLSXP: {
      const int* p = (TYPEOF(x) == INTSXP) ? INTEGER(x) : LOGICAL(x);
      for (R_xlen_t i = 1; i < n; ++i)
        if (!brk[i] && p[po[i] - 1] != p[po[i - 1] - 1]) brk[i] = 1;
      break;
    }
    case REALSXP: {
      const double* p = REAL(x);
      for (R_xlen_t i = 1; i < n; ++i) {
        if (brk[i]) continue;
        const double a = p[po[i] - 1], b = p[po[i - 1] - 1];
        if (!(a == b || (ISNAN(a) && ISNAN(b)))) brk[i] = 1;
      }
      break;
    }
    case STRSXP: {
      for (R_xlen_t i = 1; i < n; ++i) {
        if (brk[i]) continue;
        SEXP a = STRING_ELT(x, po[i] - 1), b = STRING_ELT(x, po[i - 1] - 1);
        if (a == b) continue;
        if (a == NA_STRING || b == NA_STRING || std::strcmp(CHAR(a), CHAR(b)) != 0) brk[i] = 1;
      }
      break;
    }
    default:
      stop("group_starts_cpp: unsupported key type");
    }
  }
  R_xlen_t ng = 0;
  for (R_xlen_t i = 0; i < n; ++i) ng += brk[i];
  IntegerVector out(ng);
  for (R_xlen_t i = 0, k = 0; i < n; ++i) if (brk[i]) out[k++] = (int)(i + 1);
  return out;
}

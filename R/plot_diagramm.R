#' plot a Kin diagram (network)

#' @description Given estimation of kin counts from `kins` function, draw a diagramm.
#' @param kins_total vector. Colum names must be as `kins` output.
#' @param ego_age integer.
#' @param rounding numeric. Estimation could have a lot of decimals. Rounding will make looks more clear the diagramm.
#' @return A plot
#' @importFrom DiagrammeR mermaid
#' @export


plot_diagram <- function(kins_total, ego_age = NULL, rounding = 3){
  # https://cran.r-project.org/web/packages/DiagrammeR/vignettes/graphviz-mermaid.html
  # https://color.hailpixel.com/#D9E9BE,BF62CB,94C2DB,79D297,CDA76A,C8695B
  kins_total <- round(kins_total,3)

  mermaid(
  paste0("graph TD

  GGM(ggm: <br>",                       kins_total["ggm"]  ,")
  GGM ==> GM(gm: <br>",                 kins_total["gm"]  ,")
  GM  --> AOM(oa: <br>",                kins_total["oa"]  ,")
  GM  ==> M(m: <br>",                   kins_total["m"]  ,")
  GM  --> AYM(ya: <br>",                kins_total["ya"]  ,")
  AOM  --> CAOM(coa: <br>",             kins_total["coa"]  ,")
  M   --> OS(os: <br>",                 kins_total["os"]  ,")
  M   ==> E((Ego))
  M   --> YS(ys: <br>",                 kins_total["ys"]  ,")
  AYM  --> CAYM(cya: <br>",             kins_total["cya"]  ,")
  OS   --> NOS(nos: <br>",              kins_total["nos"] ,")
  E   ==> D(d: <br>",                   kins_total["d"]  ,")
  YS   --> NYS(nys: <br>",              kins_total["nys"]  ,")
  D   ==> GD(gd: <br>",                 kins_total["gd"]  ,")
  style GGM fill:#D9E9BE, stroke:#333, stroke-width:2px;
  style GM fill:#BF62CB, stroke:#333, stroke-width:2px, text-align: center;
  style M fill:#94C2DB, stroke:#333, stroke-width:2px, text-align: center
  style D fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style YS fill:#79D297, stroke:#333, stroke-width:2px, text-align: center
  style OS fill:#79D297, stroke:#333, stroke-width:2px, text-align: center
  style CAOM fill:#79D297, stroke:#333, stroke-width:2px, text-align: center
  style AYM fill:#94C2DB, stroke:#333, stroke-width:2px, text-align: center
  style AOM fill:#94C2DB, stroke:#333, stroke-width:2px, text-align: center
  style CAYM fill:#79D297, stroke:#333, stroke-width:2px, text-align: center
  style NOS fill:#CDA76A, stroke:#333, stroke-width:2px, text-align: center
  style NYS fill:#CDA76A, stroke:#333, stroke-width:2px, text-align: center
  style E fill:#FFF, stroke:#333, stroke-width:4px, text-align: center
  style D fill:#CDA76A, stroke:#333, stroke-width:2px, text-align: center
  style GD fill:#C8695B, stroke:#333, stroke-width:2px, text-align: center"))

}



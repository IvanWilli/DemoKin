# Proposal: preserve period deaths when only one period is requested.
#
# This script is intentionally outside R/ so it is not loaded as package code.
# It shows the output-section change I would make in kin_time_variant().
#
# Problem:
#   Deaths experienced in period t are held in kin_all[[t + 1]] and then
#   shifted back by year - year_step and age_focal - 1. If kin_all is filtered
#   to only requested years before reshaping, deaths for the requested period
#   are dropped unless the following period is also requested.
#
# Main idea:
#   1. Compute the requested year-age combinations as now.
#   2. Extract requested years plus the following projection year.
#   3. Reshape living/dead values.
#   4. Relocate deaths back.
#   5. Filter final rows to the originally requested year-age combinations.

reshape_kin_time_variant_output_proposal <- function(kin_all,
                                                     out_selected,
                                                     output_kin = NULL,
                                                     age,
                                                     year_step) {
  possible_kin <- c(
    "d", "gd", "ggd", "m", "gm", "ggm", "os", "ys",
    "nos", "nys", "oa", "ya", "coa", "cya"
  )

  if (is.null(output_kin)) {
    selected_kin_position <- seq_along(possible_kin)
  } else {
    selected_kin_position <- which(possible_kin %in% output_kin)
  }

  # Include the next projection year because deaths for year t are stored in
  # the model state labelled t + year_step before relocation.
  years_requested <- unique(out_selected$year)
  years_needed <- unique(c(years_requested, years_requested + year_step))

  kin_list <- kin_all %>%
    purrr::keep(names(.) %in% as.character(years_needed)) %>%
    purrr::map(~ .[selected_kin_position])

  ages <- length(age)

  kin <- lapply(names(kin_list), FUN = function(Y) {
    X <- kin_list[[Y]]
    X <- purrr::map2(X, names(X), function(x, y) {
      x <- as.data.frame(x)
      x$year <- Y
      x$kin <- y
      x$age_kin <- rep(age, 2)
      x$alive <- c(rep("living", ages), rep("dead", ages))
      x
    }) %>%
      data.table::rbindlist() %>%
      stats::setNames(c(as.character(age), "year", "kin", "age_kin", "alive")) %>%
      data.table::melt(
        id.vars = c("year", "kin", "age_kin", "alive"),
        variable.name = "age_focal",
        value.name = "count"
      )

    X$age_focal <- as.integer(as.character(X$age_focal))
    X$year <- as.integer(X$year)
    X$cohort <- X$year - X$age_focal

    data.table::dcast(
      X,
      year + kin + age_kin + age_focal + cohort ~ alive,
      value.var = "count"
    )
  }) %>%
    data.table::rbindlist()

  living_part <- kin[, setdiff(names(kin), "dead")]
  death_part <- transform(
    kin[, setdiff(names(kin), "living")],
    year = year - year_step,
    age_focal = age_focal - 1
  )
  death_part$cohort <- death_part$year - death_part$age_focal

  kin <- merge(living_part, death_part, all.x = TRUE)

  # Now restrict to exactly the originally requested period/cohort-age rows.
  out_selected_filter <- transform(out_selected, age_focal = age)
  out_selected_filter$age <- NULL

  merge(
    out_selected_filter,
    kin,
    by = c("year", "age_focal"),
    all.x = TRUE
  )
}


# Minimal regression check to run after wiring the proposal into kin_time_variant():
#
# devtools::load_all()
# one <- kin_time_variant(
#   p = swe_px[, 1:3],
#   f = swe_asfr[, 1:3],
#   output_period = 1900,
#   output_kin = "m"
# )
# two <- kin_time_variant(
#   p = swe_px[, 1:3],
#   f = swe_asfr[, 1:3],
#   output_period = c(1900, 1901),
#   output_kin = "m"
# )
# stopifnot(any(one$dead > 0, na.rm = TRUE))
# stopifnot(isTRUE(all.equal(
#   one[one$year == 1900, c("kin", "age_kin", "age_focal", "living", "dead")],
#   two[two$year == 1900, c("kin", "age_kin", "age_focal", "living", "dead")]
# )))

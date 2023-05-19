#' rename kin

#' @description Add kin labels depending the sex
#' @details See table `demokin_codes` to know label options.
#' @param df data.frame. A data frame with variable `kin` with `DemoKin` codes to be labelled.
#' @param sex character. "f" for female, "m" for male or "2sex" for both sex naming.
#' @export
rename_kin <- function(df, sex = "f"){
  if(sex == "f") demokin_codes_sex <- DemoKin::demokin_codes[,c("DemoKin", "Labels_female")]
  if(sex == "m") demokin_codes_sex <- DemoKin::demokin_codes[,c("DemoKin", "Labels_male")]
  if(sex == "2sex") demokin_codes_sex <- DemoKin::demokin_codes[,c("DemoKin", "Labels_2sex")]
  colnames(demokin_codes_sex) <- c("kin", "kin_label")
  df %>%
    dplyr::left_join(demokin_codes_sex)
}

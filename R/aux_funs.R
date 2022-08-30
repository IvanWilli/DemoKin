
#' print kin codes
#' @description Print kin codes and labels
#' @export
demokin_codes <- function(){
  codes <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
  caswell_codes <- c("t", "v", "a", "b", "c", "h", "g", "d", "p", "q", "r", "s", "m", "n")
  labels <- c("Cousins from older aunt", "Cousins from younger aunt", "Daughter", "Grand-daughter", "Great-grand-daughter", "Great-grandmother", "Grandmother", "Mother", "Nieces from older sister", "Nieces from younger sister", "Aunt older than mother", "Aunt younger than mother", "Older sister", "Younger sister")
  data.frame(DemoKin = codes, Caswell = caswell_codes, Label = labels, row.names = NULL)
}

#' rename kin
#' @description Rename kin labels depending consolidate some types
#' @export
rename_kin <- function(df, consolidate = FALSE){

  if(!consolidate){

    relatives <- c("Cousins from older aunt", "Cousins from younger aunt", "Daughter", "Grand-daughter", "Great-grand-daughter", "Great-grandmother", "Grandmother", "Mother", "Nieces from older sister", "Nieces from younger sister", "Aunt older than mother", "Aunt younger than mother", "Older sister", "Younger sister")
    names(relatives) <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")

  } else if(consolidate){

    # Combine kin types irrespective of whether they come from older
    # or younger sibling lines
    consolidate <- c("c", "c", "d", "gd", "ggd", "ggm", "gm", "m", "n", "n", "a", "a", "s", "s")
    names(consolidate) <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")

    # Rename kin types from codes to actual words
    relatives <- c("Cousins", "Daughter", "Grand-daughter", "Great-grand-daughter", "Great-grandmother", "Grandmother", "Mother", "Nieces", "Aunt", "Sister")
    names(relatives) <-  unique(consolidate)

    df <-
      df %>%
      dplyr::mutate(kin = consolidate[kin]) %>%
      dplyr::group_by(age_focal, kin) %>%
      dplyr::summarise(count = sum(count)) %>%
      dplyr::ungroup()

  }
  df$kin <- relatives[df$kin]
  df
}

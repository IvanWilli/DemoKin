
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
rename_kin <- function(df, consolidate_column = "no"){

  stopifnot("Argument 'consolidate_column' should be 'no' or a valid column name" = consolidate_column %in% c("no", colnames(df)))

  if(consolidate_column == "no"){

    relatives <- c("Cousins from older aunt", "Cousins from younger aunt", "Daughter", "Grand-daughter", "Great-grand-daughter", "Great-grandmother", "Grandmother", "Mother", "Nieces from older sister", "Nieces from younger sister", "Aunt older than mother", "Aunt younger than mother", "Older sister", "Younger sister")
    names(relatives) <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")

  } else {

    # Combine kin types irrespective of whether they come from older
    # or younger sibling lines
    consolidate_vec <- c("c", "c", "d", "gd", "ggd", "ggm", "gm", "m", "n", "n", "a", "a", "s", "s")
    names(consolidate_vec) <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")

    # Rename kin types from codes to actual words
    relatives <- c("Cousins", "Daughter", "Grand-daughter", "Great-grand-daughter", "Great-grandmother", "Grandmother", "Mother", "Nieces", "Aunt", "Sister")
    names(relatives) <-  unique(consolidate_vec)

    df <- as.data.frame(df)
    df$count <- df[ , consolidate_column]

    df <-
      df %>%
      dplyr::mutate(kin = consolidate_vec[kin]) %>%
      dplyr::group_by(age_focal, kin) %>%
      dplyr::summarise(count = sum(count)) %>%
      dplyr::ungroup()


  }
  df$kin <- relatives[df$kin]
  df
}

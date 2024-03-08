#' plot a Kin diagram (network)

#' @description Draws a Keyfitz-style kinship diagram given a kinship object created by the `kin` function. Displays expected kin counts for a Focal aged 'a'.
#' @param kin_total data.frame. values in column `kin` define the relative type - see `demokin_codes()`. Values in column `count` are the expected number of relatives.
#' @param rounding numeric. Number of decimals to show in diagram.
#' @return A Keyfitz-style kinship plot.
#' @export

plot_diagram <-
  function (kin_total, rounding = 3) {
    rels <- c("ggd", "gd", "d", "Focal", "m", "gm", "ggm", "oa", "coa", "os", "nos", "ya", "cya", "ys", "nys")
    # check all types are in
    if(!any(unique(kin_total$kin) %in% rels) | any(c("s", "c", "a", "n") %in% unique(kin_total$kin))) stop("You need all specific types. If some are missed or grouped, for example old and younger sisters in 's', this will fail.")
    vertices <- data.frame(
      nodes = rels
      , x = c(1, 1, 1, 1, 1, 1, 1, 0, -1, 0, -1, 2, 3, 2, 3)
      , y = c(0, 1, 2, 3, 4, 5, 6, 4, 3, 3, 2, 4, 3, 3, 2)
    )
    d <- data.frame(from = c("ggd", "gd", "d", "Focal", "m",
                             "gm", "gm", "oa", "m", "os", "gm", "ya", "m", "ys"),
                    to = c("gd", "d", "Focal", "m", "gm", "ggm", "oa", "coa",
                           "os", "nos", "ya", "cya", "ys", "nys"))
    lookup <- c(with(kin_total, paste0(kin, " \n", round(count, rounding))), "Focal")
    names(lookup) <- c(kin_total$kin, "Focal")
    vertices$nodes <- lookup[vertices$nodes]
    d$from <- lookup[d$from]
    d$to <- lookup[d$to]
    # to show full relative names
    relatives <- c("Cousins from older aunt", "Cousins from younger aunt",
                   "Daughter", "Grand-daughter", "Great-grand-daughter",
                   "Great-grandmother", "Grandmother", "Mother", "Nieces from older sister",
                   "Nieces from younger sister", "Aunt older than mother",
                   "Aunt younger than mother", "Older sister", "Younger sister", "")
    names(relatives) <- c("coa", "cya", "d", "gd", "ggd",
                          "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os",
                          "ys", "Focal")
    labs <- relatives[rels]
    # Plot
    b <- igraph::graph_from_data_frame(vertices = vertices, d= d, directed = FALSE)
    b_auto_layout <- igraph::layout.auto(b)
    b_auto_layout_scaled <- igraph::norm_coords(b_auto_layout, ymin=-1, ymax=1, xmin=-1, xmax=1)
    plot(
      b
      , vertex.size = 70
      , curved = 1
      , vertex.color = "#FFF1E2"
      , vertex.shape = "circle"
      , vertex.label.cex = 0.8
      , vertex.label.color = "black"
      , edge.width = 2
      , layout = b_auto_layout_scaled * 3
      , rescale = FALSE
      , xlim = c(-3.3,3.3)
      , ylim = c(-3.1,3.1)
    )
    # Add relative names
    # Thanks to Egor Kotov for this tip!
    plot(
      b
      , vertex.size = 70
      , curved = 1
      , vertex.color = NA
      , vertex.shape = "none"
      , vertex.label = labs
      , vertex.label.dist = -6.5
      , vertex.label.cex = 0.8
      , vertex.label.color = "black"
      , vertex.label.degree = -pi/2
      , edge.width = 2
      , edge.color = NA
      , layout = b_auto_layout_scaled * 3
      , rescale = FALSE
      , xlim = c(-3.3,3.3)
      , ylim = c(-3.1,3.1)
      , add = T
    )
  }

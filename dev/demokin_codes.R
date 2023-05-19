#' demokin kin codes

codes <- c("coa", "cya", 'c', "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "n", "oa", "ya", "a","os", "ys", "s")
caswell_codes <- c("t", "v", NA, "a", "b", "c", "h", "g", "d", "p", "q", NA, "r", "s", NA,"m", "n",NA)
labels_female <- c("Cousins from older aunts", "Cousins from younger aunts", "Cousins",
                   "Daughters", "Grand-daughters",
                   "Great-grand-daughters", "Great-grandmothers", "Grandmothers", "Mother",
                   "Nieces from older sisters", "Nieces from younger sisters",  "Nieces",
                   "Aunts older than mother", "Aunts younger than mother",  "Aunts",
                   "Older sisters", "Younger sisters", "Sisters")
labels_male <- c("Cousins from older uncles", "Cousins from younger uncles", "Cousins",
                 "Brothers", "Grand-sons",
                 "Great-grand-sons", "Great-grandfathers", "Grandfathers", "Father",
                 "Nephews from older brothers", "Nephews from younger brothers", "Nephews",
                 "Uncles older than fathers", "Uncles younger than father", "Uncles",
                 "Older brothers", "Younger brothers", "Brothers")
labels_2sex <- c("Cousins from older aunts/uncles", "Cousins from younger aunts/uncles", "Cousins",
                 "Siblings", "Grand-childrens",
                 "Great-grand-childrens", "Great-grandfparents", "Grandparents", "Parents",
                 "Niblings from older siblings", "Niblings from younger siblings", "Niblings",
                 "Aunts/Uncles older than parents", "Aunts/Uncles younger than parents", "Aunts/Uncles",
                 "Older siblings", "Younger siblings", "Siblings")
demokin_codes <- data.frame(DemoKin = codes, Caswell = caswell_codes,
           Labels_female = labels_female,
           Labels_male   = labels_male,
           Labels_2sex   = labels_2sex,
           row.names = NULL)
# save(demokin_codes, file = "data/demokin_codes.rda")

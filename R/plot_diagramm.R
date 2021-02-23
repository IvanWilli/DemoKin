# graph -------------------------------------------------------------------
plot_diagramm <- function(kins_total){
  # https://cran.r-project.org/web/packages/DiagrammeR/vignettes/graphviz-mermaid.html

  mermaid(
    paste0("graph TB
  GGM(ggmother: <br>",              round(kins_total["ggm"],5)  ,")
  GGM ==> GM(gmother: <br>",        round(kins_total["gm"],5)  ,")
  GM  --> AOM(aunt older mother: <br>", .5  ,")
  GM  ==> M(mother: <br>",          round(kins_total["m"],5)  ,")
  GM  --> AYM(aunt younger mother: <br>", .5  ,")
  AOM  --> CAOM(cousin from aunt older mother: <br>", .5  ,")
  M   --> OS(older sister: <br>",   round(kins_total["os"],5)  ,")
  M   ==> E((ego: <br>", 1  ,"))
  M   --> YS(younger sis: <br>",    round(kins_total["ys"],5)  ,")
  AYM  --> CAYM(cousin from aunt younger mother: <br>", .5  ,")
  OS   --> NOS(niece from older sister: <br>", .5  ,")
  E   ==> D(daughter: <br>",        round(kins_total["d"],5)  ,")
  YS   --> NYS(niece from younger sister: <br>", .5  ,")
  D   ==> GD(ggdaughter: <br>",     round(kins_total["gd"],5)  ,")
  style GGM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style GM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style M fill:#ffe8e2, stroke:#333, stroke-width:2px
  style D fill:#ffe8e2, stroke:#333, stroke-width:2px
  style YS fill:#ffe8e2, stroke:#333, stroke-width:2px
  style OS fill:#ffe8e2, stroke:#333, stroke-width:2px
  style CAOM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style AYM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style AOM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style CAYM fill:#ffe8e2, stroke:#333, stroke-width:2px
  style NOS fill:#ffe8e2, stroke:#333, stroke-width:2px
  style NYS fill:#ffe8e2, stroke:#333, stroke-width:2px
  style E fill:#FFF, stroke:#333, stroke-width:4px
  style D fill:#ffe8e2, stroke:#333, stroke-width:2px
  style GD fill:#ffe8e2, stroke:#333, stroke-width:2px"))
}



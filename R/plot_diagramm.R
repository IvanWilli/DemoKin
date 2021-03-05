# graph -------------------------------------------------------------------
plot_diagram <- function(kins_total, ego_age = NULL, rounding = 3){
  # https://cran.r-project.org/web/packages/DiagrammeR/vignettes/graphviz-mermaid.html

  kins_total <- round(kins_total,3)

  mermaid(
  paste0("graph TB

  title[<u>Expected kins for an ego aged ", ego_age ," </u>]
  title-->GGM
  style title fill:#FFF,stroke:#FFF,font-size:150%;
  linkStyle 0 stroke:#FFF,stroke-width:0;

  GGM(Great-grandgmother: <br>",                      round(kins_total["ggm"],5)  ,")
  GGM ==> GM(Grandmother: <br>",                round(kins_total["gm"],5)  ,")
  GM  --> AOM(Aunts older than Mother: <br>",     round(kins_total["oa"],5)  ,")
  GM  ==> M(Mother: <br>",                    round(kins_total["m"],5)  ,")
  GM  --> AYM(Aunts younger than mother: <br>",   round(kins_total["ya"],5)  ,")
  AOM  --> CAOM(Cousins: <br>",               round(kins_total["coa"],5)  ,")
  M   --> OS(Older sister: <br>",             round(kins_total["os"],5)  ,")
  M   ==> E((Ego: <br>",                        1  ,"))
  M   --> YS(Younger sisters: <br>",            round(kins_total["ys"],5)  ,")
  AYM  --> CAYM(Cousins: <br>",               round(kins_total["cya"],5)  ,")
  OS   --> NOS(Nieces through older sister: <br>",   round(kins_total["nos"],5)  ,")
  E   ==> D(Daughters: <br>",                 round(kins_total["d"],5)  ,")
  YS   --> NYS(Nieces through younger sister: <br>", round(kins_total["nys"],5)  ,")
  D   ==> GD(Granddaughters: <br>",             round(kins_total["gd"],5)  ,")
  style GGM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center;
  style GM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center;
  style M fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style D fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style YS fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style OS fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style CAOM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style AYM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style AOM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style CAYM fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style NOS fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style NYS fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style E fill:#FFF, stroke:#333, stroke-width:4px, text-align: center
  style D fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center
  style GD fill:#dddbdb, stroke:#333, stroke-width:2px, text-align: center"))
}



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

  kins_total <- round(kins_total,3)

  mermaid(
  paste0("graph TB

  title[<u>Expected kins for an ego aged ", ego_age ," </u>]
  title-->GGM
  style title fill:#FFF,stroke:#FFF,font-size:150%;
  linkStyle 0 stroke:#FFF,stroke-width:0;

  GGM(ggm: <br>",                      round(kins_total["ggm"],5)  ,")
  GGM ==> GM(gm: <br>",                round(kins_total["gm"],5)  ,")
  GM  --> AOM(oa: <br>",     round(kins_total["oa"],5)  ,")
  GM  ==> M(m: <br>",                    round(kins_total["m"],5)  ,")
  GM  --> AYM(ya: <br>",   round(kins_total["ya"],5)  ,")
  AOM  --> CAOM(coa: <br>",               round(kins_total["coa"],5)  ,")
  M   --> OS(os: <br>",             round(kins_total["os"],5)  ,")
  M   ==> E((Ego: <br>",                        1  ,"))
  M   --> YS(ys: <br>",            round(kins_total["ys"],5)  ,")
  AYM  --> CAYM(cya: <br>",               round(kins_total["cya"],5)  ,")
  OS   --> NOS(nos: <br>",   round(kins_total["nos"],5)  ,")
  E   ==> D(d: <br>",                 round(kins_total["d"],5)  ,")
  YS   --> NYS(nys: <br>", round(kins_total["nys"],5)  ,")
  D   ==> GD(gd: <br>",             round(kins_total["gd"],5)  ,")
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



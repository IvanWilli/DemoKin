#' plot a Kin diagram (network)

#' @description Given estimation of kin counts from `kins` function, draw a diagramm.
#' @param kins_total vector. Colum names must be as `kins` output.
#' @param rounding numeric. Estimation could have a lot of decimals. Rounding will make looks more clear the diagramm.
#' @return A plot
#' @importFrom DiagrammeR mermaid
#' @export


plot_diagram <- function(kins_total, rounding = 3){
  # https://cran.r-project.org/web/packages/DiagrammeR/vignettes/graphviz-mermaid.html
  # https://color.hailpixel.com/#D9E9BE,BF62CB,94C2DB,79D297,CDA76A,C8695B

  kins_total <- kins_total %>% mutate(total = round(total,digits = rounding))

  mermaid(
  paste0("graph TD

  GGM(ggm: <br>",                       kins_total$total[kins_total$kin=="ggm"]  ,")
  GGM ==> GM(gm: <br>",                 kins_total$total[kins_total$kin=="gm"]  ,")
  GM  --> AOM(oa: <br>",                kins_total$total[kins_total$kin=="oa"]  ,")
  GM  ==> M(m: <br>",                   kins_total$total[kins_total$kin=="m"]  ,")
  GM  --> AYM(ya: <br>",                kins_total$total[kins_total$kin=="ya"]  ,")
  AOM  --> CAOM(coa: <br>",             kins_total$total[kins_total$kin=="coa"]  ,")
  M   --> OS(os: <br>",                 kins_total$total[kins_total$kin=="os"]  ,")
  M   ==> E((Ego))
  M   --> YS(ys: <br>",                 kins_total$total[kins_total$kin=="ys"]  ,")
  AYM  --> CAYM(cya: <br>",             kins_total$total[kins_total$kin=="cya"]  ,")
  OS   --> NOS(nos: <br>",              kins_total$total[kins_total$kin=="nos"] ,")
  E   ==> D(d: <br>",                   kins_total$total[kins_total$kin=="d"]  ,")
  YS   --> NYS(nys: <br>",              kins_total$total[kins_total$kin=="nys"]  ,")
  D   ==> GD(gd: <br>",                 kins_total$total[kins_total$kin=="gd"]  ,")
  style GGM fill:#a1f590, stroke:#333, stroke-width:2px;
  style GM  fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center;
  style M   fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center
  style D   fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center
  style YS  fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center
  style OS  fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center
  style CAOM fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style AYM fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style AOM fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style CAYM fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style NOS fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style NYS fill:#f1f0f5, stroke:#333, stroke-width:2px, text-align: center
  style E   fill:#FFF, stroke:#333, stroke-width:4px, text-align: center
  style D   fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center
  style GD  fill:#a1f590, stroke:#333, stroke-width:2px, text-align: center"))

}



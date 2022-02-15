#' plot a Kin diagram (network)

#' @description Given estimation of kin counts from `kins` function, draw a diagramm.
#' @param kins_total vector. Colum names must be as `kins` output.
#' @param rounding numeric. Estimation could have a lot of decimals. Rounding will make looks more clear the diagramm.
#' @return A plot
#' @importFrom DiagrammeR mermaid
#' @export


plot_diagram <- function(kin_total, rounding = 3){
  # https://cran.r-project.org/web/packages/DiagrammeR/vignettes/graphviz-mermaid.html
  # https://color.hailpixel.com/#D9E9BE,BF62CB,94C2DB,79D297,CDA76A,C8695B

  kin_total <- kin_total %>% mutate(count = round(count,digits = rounding))

  mermaid(
  paste0("graph TD

  GGM(ggm: <br>",                       kin_total$count[kin_total$kin=="ggm"]  ,")
  GGM ==> GM(gm: <br>",                 kin_total$count[kin_total$kin=="gm"]  ,")
  GM  --> AOM(oa: <br>",                kin_total$count[kin_total$kin=="oa"]  ,")
  GM  ==> M(m: <br>",                   kin_total$count[kin_total$kin=="m"]  ,")
  GM  --> AYM(ya: <br>",                kin_total$count[kin_total$kin=="ya"]  ,")
  AOM  --> CAOM(coa: <br>",             kin_total$count[kin_total$kin=="coa"]  ,")
  M   --> OS(os: <br>",                 kin_total$count[kin_total$kin=="os"]  ,")
  M   ==> E((Ego))
  M   --> YS(ys: <br>",                 kin_total$count[kin_total$kin=="ys"]  ,")
  AYM  --> CAYM(cya: <br>",             kin_total$count[kin_total$kin=="cya"]  ,")
  OS   --> NOS(nos: <br>",              kin_total$count[kin_total$kin=="nos"] ,")
  E   ==> D(d: <br>",                   kin_total$count[kin_total$kin=="d"]  ,")
  YS   --> NYS(nys: <br>",              kin_total$count[kin_total$kin=="nys"]  ,")
  D   ==> GD(gd: <br>",                 kin_total$count[kin_total$kin=="gd"]  ,")
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



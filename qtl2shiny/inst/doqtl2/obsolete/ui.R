## https://gist.github.com/wch/4211337
## dynamic input fields
## http://shiny.stat.ubc.ca/r-graph-catalog/

library(shiny)
#library(shinyAce)

## Function from Joe Cheng
## https://gist.github.com/jcheng5/5913297
helpPopup <- function(title, content,
                      placement = c('right', 'top', 'left', 'bottom'),
                      trigger = c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover()})"),
        tags$style(type = "text/css", ".popover{max-width:500px; position: fixed;}")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-link",
      `data-toggle` = "popover", `data-html` = "true",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok = TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok = TRUE)[1],
      "More..."
    )
  )
}

## Function adapted from Jenny Bryan
## https://github.com/jennybc/r-graph-catalog
shinyUI(navbarPage("DO QTL2",

  tabPanel("Phenos",
    shinyPhenosUI("phenos"),
    shinyWindowUI("window")
  ),

  tabPanel("Peaks",
           shinyPeaksUI("shinypeaks")
  ),

  tabPanel("Tables",
    h4("analyses"),
    tableOutput("analyses_tbl"),
    h4("peaks"),
    tableOutput("peaks_tbl")
  ),

  tabPanel("Plot",
    selectInput("plotType", "Plot Type",
                c("density", "scans", "snps")),
    conditionalPanel(condition="input.plotType == 'density'",
                    plotOutput("distPlot")),
    conditionalPanel(condition="input.plotType == 'scans'",
                     shinyScan1UI("genome_scan")),
    conditionalPanel(condition="input.plotType == 'snps'",
                     shinyScan1UI("snp_scan"))
),

  tabPanel("One SNP",
           actionButton("snpButton", "markdown"),
           plotOutput("heatmap"),
           uiOutput("pheno_string"),
           ## This does not quite work.
           includeMarkdown("html.md")
  ),

  tabPanel("More",
    includeMarkdown("attieDO.md"),
    includeMarkdown("about.md")
  ),
  tags$div(id = "popup",
           helpPopup(strong("Additional Information"),
                     includeMarkdown("about-extended.md"),
                     placement = "right", trigger = "click"))
))


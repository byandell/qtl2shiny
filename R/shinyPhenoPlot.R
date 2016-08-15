#' Shiny Phenotype Plot module
#'
#' Shiny module to plot phenotypes.
#'
#' @param input,output,session standard shiny arguments
#' @param raw_phe_df,phe_df,cov_mx reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return 2-element vector of scan window
#'
#' @export
shinyPhenoPlot <- function(input, output, session,
                       raw_phe_df, phe_df, cov_mx) {
  ## Scatter plot or density
  mysum <- function(phe) {
    tmpfn <- function(x) {
      out <- summary(x)
      out["NA's"] <- sum(is.na(x))
      out
    }
    t(sapply(phe, tmpfn))
  }
  output$phe_sum <- renderTable({
    mysum(phe_df())
  })
  output$raw_phe_sum <- renderTable({
    mysum(raw_phe_df())
  })
  myplot <- function(phe, cov) {
    phename <- names(phe)
    if(length(phename) > 10) {
      cat(file=stderr(), "\nOnly first 10 phenotypes used\n")
      phename <- phename[seq_len(10)]
      phe <- phe[,phename]
    }
    if("sex" %in% dimnames(cov)[[2]]) {
      ## Assume sex in covar. Ignore actual covariates for analyses.
      insex <- data.frame(phe,cov) %>%
        mutate(sex=c("female","male")[1+sex])

      if(length(phename) == 1) {
        ggplot(insex, aes_string(phename, col="sex")) +
          geom_density() + geom_rug()
      } else {
        ggscatmat(insex, seq_along(phe), color="sex")
      }
    } else {
      if(length(phename) == 1) {
        ggplot(phe, aes_string(phename)) +
          geom_density() + geom_rug()
      } else {
        ggscatmat(phe, seq_along(phe))
      }
    }
  }
  output$phePlot <- renderPlot({
    myplot(phe_df(), cov_mx())
  })
  output$rawphePlot <- renderPlot({
    myplot(raw_phe_df(), cov_mx())
  })
}

#' UI for shinyPhenoPlot Shiny Module
#'
#' UI for phenotype plots and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyPhenoPlot
#' @export
shinyPhenoPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("raw pheno",
               plotOutput(ns("rawphePlot")),
               tableOutput(ns("raw_phe_sum"))),
      tabPanel("transformed",
               plotOutput(ns("phePlot")),
               tableOutput(ns("phe_sum")))
    )
  )
}

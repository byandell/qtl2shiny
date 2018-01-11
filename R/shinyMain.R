#' Shiny app for DOQTL selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param projects_info reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom gdata humanReadable
#' @importFrom dplyr distinct filter mutate one_of select 
#' @importFrom tidyr unite
#' @importFrom qtl2pattern pheno_trans
#' @importFrom shiny callModule NS reactive req 
#'   radioButtons selectInput
#'   dataTableOutput textOutput uiOutput
#'   renderDataTable renderText renderUI
#'   observeEvent
#'   strong tagList
shinyMain <- function(input, output, session, projects_info) {
  ns <- session$ns
  
  ## Data Setup
  project_info <- reactive({ 
    shiny::req(set_par()$project_info) 
  })
  pheno_type <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "pheno_type")
  })
  peaks <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "peaks")
  })
  analyses_tbl <- shiny::reactive({
    shiny::req(project_info())
    ## The analyses_tbl should only have one row per pheno.
    read_project_rds(project_info(), "analyses")
  })
  pheno_data <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "pheno_data")
  })
  covar <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "covar")
  })
  
  pmap_obj <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "pmap")
  })
  kinship <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "kinship")
  })
  
  set_par <- shiny::callModule(shinySetup, "setup", 
                               pheno_type, peaks, 
                               pmap_obj, analyses_tbl, 
                               cov_df, pheno_data, projects_info)
  
  ## Continue with Plots and Analysis.
  
  ## Phenotypes and Covariates.
  analyses_df <- shiny::reactive({
    phename <- shiny::req(set_par()$phe_par$pheno_names)
    dplyr::filter(analyses_tbl(), pheno %in% phename)
  })
  phe_mx <- shiny::reactive({
    shiny::req(analyses_df())
    pheno_read(pheno_data(), analyses_df())
  })
  cov_df <- shiny::reactive({
    DOread::get_covar(covar(), analyses_df())
  })
  
  ## Set up shiny::reactives for scan1 module.
  K_chr <- shiny::reactive({
    kinship()[set_par()$win_par$chr_id]
  })
  
  ## Allele names.
  allele_info <- shiny::reactive({
    shiny::req(project_info())
    read_project_rds(project_info(), "allele_info")
  })
  
  ## Haplotype Analysis.
  shiny::callModule(shinyHaplo, "hap_scan", 
                    set_par()$win_par, pmap_obj, 
                    phe_mx, cov_df, K_chr, analyses_df, 
                    covar, pheno_data, analyses_tbl, peaks,
                    project_info, allele_info)
  
  ## Diplotype Analysis.
  shiny::callModule(shinyDiplo, "dip_scan",
                    set_par()$win_par, 
                    phe_mx, cov_df, K_chr, analyses_df,
                    project_info, allele_info)
}

#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyMain
#' @export
shinyMainInput <- function(id) {
  ns <- shiny::NS(id)
  shinySetupInput(ns("setup"))
}
#' @rdname shinyMain
#' @export
shinyMainUI <- function(id) {
  ns <- shiny::NS(id)
  shinySetupUI(ns("setup"))
}
#' @rdname shinyMain
#' @export
shinyMainOutput <- function(id) {
  ns <- shiny::NS(id)
  shinyHaploUI(ns("hap_scan"))
}
#' @rdname shinyMain
#' @export
shinyMainOutput2 <- function(id) {
  ns <- shiny::NS(id)
  shinyDiploUI(ns("dip_scan"))
}

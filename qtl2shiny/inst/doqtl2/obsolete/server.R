suppressPackageStartupMessages({
  library(doqtl2)
  library(shiny)
  library(GGally) #ggscatmat
  library(markdown) #includeMarkdown()
})

dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")
pmap <- readRDS(file.path(datapath, "pmap.rds"))
peaks <- readRDS(file.path(datapath, "peaks.rds"))
covar <- readRDS(file.path(datapath, "covar.rds"))
peak_info <- peaks$output
analyses_tbl <- readRDS(file.path(datapath, "analyses.rds")) %>%
  filter(output %in% peak_info)
peak_info <- peaks$pheno
pheno_data <- read_pheno_tbl(analyses_tbl, datapath)
pheno_data <- pheno_data %>%
  select(which(names(pheno_data) %in% peak_info))
rm(peak_info)

pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))
K <- readRDS(file.path(datapath, "kinship.rds"))

shinyServer(function(input, output, session) {
  ##############################################
  ## Set up phenotype choice.

  pmap_obj <- reactive({pmap})

  win_par <- callModule(shinyWindow, "window", pmap_obj)
  chr_id <- reactive({win_par$chr_id})
  peak_Mbp <- reactive({win_par$peak_Mbp})
  window_Mbp <- reactive({win_par$window_Mbp})

  ## Reactives for shiny modules.
  pheno_typer <- reactive({pheno_type})
  analyses_tblr <- reactive({analyses_tbl})
  peaks_tbl <- reactive({peaks})

  ## Call shiny module to look for peaks.
  ## Want output from this to automate peak selection. Later.
  callModule(shinyPeaks, "shinypeaks",
             pheno_typer, peaks_tbl, pmap_obj)

  ## Call shiny module to update phenotype list.
  pheno_re <- callModule(shinyPhenos, "phenos",
                         pheno_typer, peaks_tbl, analyses_tblr,
                         win_par)

  ##############################################
  ## Continue with Plots and Analysis.

  ## Set up peaks data frame.
  peaks_df <- reactive({
    peaks %>%
      filter(output %in% pheno_re())
  })

  ## Set up analyses data frame.
  analyses_df <- reactive({
    analyses_tbl %>%
      filter(output %in% pheno_re())
  })

  # Output the analyses table
  output$analyses_tbl <- renderTable({
    dat <- analyses_df()

    ## Collapse covariates (past winsorize column).
    covar_names <- names(dat)[-seq_len(match("winsorize",names(dat)))]
    dat %>%
      unite(covar, one_of(covar_names)) %>%
      mutate(covar = sapply(strsplit(covar,"_"),
                            function(x) paste(covar_names[as.logical(x)],
                                              collapse=",")))
  })

  # Output the peaks table
  output$peaks_tbl <- renderTable({peaks_df()})

  ## Reactive for phenotypes.
  phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno))
  })

  ## Reactive for covariates
  cov_mx <- reactive({
    get_covar(covar, analyses_df())
  })

  ## Scatter plot when button pushed
  output$distPlot <- renderPlot({
    phe <- phe_df()
    phename <- names(phe)
    if(length(phename) < 2)
      return()
    if(length(phename) > 10) {
      cat(file=stderr(), "\nOnly first 10 phenotypes used\n")
      phename <- phename[seq_len(10)]
      phe <- phe[,phename]
    }

    ## Assume sex in covar. Ignore actual covariates for analyses.
    insex <- data.frame(phe,covar) %>%
      mutate(sex=c("female","male")[1+sex])

    if(length(phename) == 1) {
      p <- ggplot(insex, aes(phename[1], col="sex")) +
        geom_density()
    } else {
      p <- ggscatmat(insex, seq_along(phe), color="sex")
    }
    p
  })

  ## Set up reactives for scan1 module.
  chr_id <- reactive({as.character(req(win_par$chr_id))})
  K_chr <- reactive({K[chr_id()]})
  plot_type <- reactive({input$plotType})

  ## Genome scan.
  probs_obj <- reactive({
    read_probs(chr_id(), datapath)
  })
  callModule(shinyScan1, "genome_scan", plot_type,
             chr_id, phe_df, cov_mx, pheno_re, analyses_df,
             probs_obj, K_chr)

  ## SNP scan.
  snpprobs_obj <- reactive({
    window_Mbp <- req(win_par$window_Mbp)
    peak_Mbp <- req(win_par$peak_Mbp)
    if(window_Mbp == 0) {
      window_Mbp <- 3
      cat(file=stderr(),
          "\nNo window_Mbp provided -- set to 3\n")
    }
    if(peak_Mbp == 0) {
      cat(file=stderr(),
          "\nNo peak_Mbp provided -- set to pre-computed peak if available\n")
      peak_Mbp <- mean((peaks_df() %>%
                          filter(chr==chr_id()) %>%
                          filter(lod==max(lod)))$pos_Mbp)
      if(is.na(peak_Mbp)) {
        cat(file=stderr(),
            paste("\nNo pre-computed significant peaks on chr", chr_id(), "\nChange Params settings\n"))
      }
    }
    snpinfo <- get_snpinfo(chr_id(), peak_Mbp, window_Mbp, pattern = "AB1NZCPW", datapath)
    genoprob_to_snpprob(probs_obj(), snpinfo)
  })
  callModule(shinyScan1, "snp_scan", plot_type,
               chr_id, phe_df, cov_mx, pheno_re, analyses_df,
               snpprobs_obj, K_chr)

  ## Experiment with d3 object. Outputs to minipage.
  output$heatmap <- renderPlot({
    out_lmm_snps <- lmm_snps_obj()
    max_pos <- max(out_lmm_snps)$pos
    rng <- max_pos + c(-.5,.5)
    x <- out_lmm_snps$map[[1]]
    wh <- which(tmp >= rng[1] & tmp <= rng[2])
    lod <- out_lmm_snps$lod[wh,]
    x <- x[wh]
    y <- seq_len(ncol(lod))
    qtlcharts::iheatmap(lod, x, y)
  })
  ## Placeholder for Rmarkdown. Does not quite work.
  output$pheno_string <- renderText({
    # Take a dependency on input$goButton
    observeEvent(input$snpButton, {
      cat("[otu](file:///Users/brianyandell/Documents/Research/attie_alan/DO/doqtl2/inst/doqtl2/otu.html)\n",
          file="html.md")
      rmarkdown::render("otu.Rmd")
      "html.md"
    })

    pheno_re()
  })
})

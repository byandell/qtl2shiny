Setup 209
  project_info callModule Project
  win_par callModule Peaks
  phe_par callModule Phenos
  observeEvent project_info
    num_phenos renderText
  observeEvent phe_par$pheno_names
    num_phenos renderText
SetupInput
  num_phenos renderText
  chr_pos renderText
SetupUI
  title renderUI
    Region
      ProjectUI
    Phenotypes
      renderText strong
  radio radioButtons
  sidebar_setup renderUI
    Region
      PeaksInput
    Phenotypes
      filter renderUI checkboxInput
      pheno_names renderUI selectInput
      PhenosUI
  pheno_group selectInput
  dataset selectInput
  
  main_setup renderUI
    Region
      PeaksOutput
    Phenotypes
      PhenosOutput

ProjectUI 43
  project renderUI selectInput
  
Peaks 169
  callModule Hotspot
  observeEvent project_info
    chr_id updateSelectInput
    window_Mbp updateNumericInput
  observeEvent chr_id
    peak_Mbp updateNumericInput
  observeEvent Hotspot
    chr_id updateSelectInput
    peak_Mbp updateNumericInput
  observeEvent chr_pos
    chr_id updateSelectInput
    peak_Mbp updateNumericInput
PeaksInput
  local checkboxInput
  columns
    chr_id renderUI selectInput
    peak_Mbp renderUI numericInput
    window_Mbp renderUI numericInput
  HotspotInput
PeaksOutput
  HotspotOutput
  
Hotspot 166
  observeEvent project_info
    chr_ct updateSelectInput
  observeEvent chr_ct
    chr_ct updateSelectInput
HotspotInput
  fluidRow
    chr_ct renderUI selectInput
    minLOD renderUI numericInput
    peak_ck renderUI checkboxInput
HotspotOutput
  peak_show renderUI
    peak_plot renderPlot
  peak_tbl renderDataTable
  
Phenos 129
  callModule PhenoPlot raw
  callModule PhenoPlot trans
PhenosUI
  radio radioButtons
PhenosOutput
  show_data renderUI
    peaks_tbl renderDataTable dataTableOutput
    raw PhenoPlotUI
    trans PhenoPlotUI
    analyses_tbl renderDataTable dataTableOutput
  pheno_lod renderDataTable
  
PhenoPlot 66
  phePlotTable renderUI
PhenoPlotUI
  phePlot renderPlot
  pheSum renderTable
  
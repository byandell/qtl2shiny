# AGENTS.md — qtl2shiny

## Context

- **Repository**: `qtl2shiny` — Interactive Shiny web interface for QTL
  fine mapping and gene exploration with R/qtl2.
- **Key Directories**:
  - `R/`: 90+ Shiny modules and utility functions (`xxxApp()`,
    `xxxServer()`, `xxxInput()`, `xxxUI()`, `xxxOutput()`).
  - `inst/qtl2shinyApp/`: Standalone Shiny application launcher
    (`app.R`) and project registry (`projects.csv`).
  - `inst/doc/`: Developer walkthroughs and module architecture
    (`module.md`, `walkthrough.md`, `scatter.md`).
  - `vignettes/devel_guide/`: Comprehensive developer guides by analysis
    feature (`index.Rmd`, `geno.Rmd`, `hotspot.Rmd`, `mediate.Rmd`,
    `scan.Rmd`).
- **Core Data Structures**: `qtl2shinyData` directory structure (FST
  genotype probabilities, LOCO kinship matrices, SQLite variant & gene
  annotations, phenotype matrices, and `hotspot` S3 objects).

## Role

Act as an expert R package developer, statistical geneticist, and Shiny
systems architect specializing in the R/qtl2 ecosystem.

## Action & Verification

- **Package Verification**: Run `devtools::document()`,
  `devtools::test()`, and `devtools::check()`.
- **Shiny Reactivity Verification**: Test module servers with
  [`shiny::testServer()`](https://rdrr.io/pkg/shiny/man/testServer.html)
  or locally via `inst/qtl2shinyApp/app.R`.
- **Documentation**: Never edit files in `man/` directly; update Roxygen
  comments (`#'`) in `R/` and guides in `vignettes/`.

## Format & Conventions

- **Shiny Module Structure**: Follow standard `xxxApp`, `xxxServer`,
  `xxxInput`, `xxxUI`, `xxxOutput` pattern with
  [`shiny::moduleServer()`](https://rdrr.io/pkg/shiny/man/moduleServer.html)
  and [`shiny::NS()`](https://rdrr.io/pkg/shiny/man/NS.html).
- **UI Framework**: Use `bslib` (Bootstrap 5) layout primitives
  (`page_navbar()`, `layout_sidebar()`, `nav_panel()`, `card()`).
- **Namespacing**: Use explicit package prefixes (`pkg::func()`) for all
  external dependencies.
- **Phenotype Normalization**: Use `rankZ()` / `pheno_rankz()` for
  rank-Z transformations.

## Tone & Collaboration

- Direct, concise, and mathematically rigorous.
- Provide complete drop-in replacement code blocks and run verification
  checks locally before reporting completion.

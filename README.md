
TimeCycle <img src="man/figures/logo.png" align="right"/>
=========================================================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
Status](https://www.r-pkg.org/badges/version/TimeCycle)](https://cran.r-project.org/package=TimeCycle)
[![R build
status](https://github.com/r-lib/pkgdown/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/pkgdown/actions)
[![R build
status](https://github.com/nesscoder/TimeCycle/workflows/R-CMD-check/badge.svg)](https://github.com/nesscoder/TimeCycle/actions)
<!-- badges: end -->

> TimeCycle is designed to …

Installation
------------

    # Install release version from CRAN
    install.packages("TimeCycle")

    # Install development version from GitHub
    devtools::install_github("nesscoder/TimeCycle")

Usage
-----

This is a basic example which shows you how to solve a common problem:
See Reference&gt;zhang2014 For additional info on dataset. Set the
replicate labels.

    library(TimeCycle)

    #set seed for reproducibility with random variables in example usage
    set.seed(1234) 

    TimeCycleResults <- TimeCycle(data = zhang2014, repLabel = rep(1,24), period = 24)
    #> 
    #>       ########################################################################################
    #>       ###      ████████ ██ ███    ███ ███████  ██████ ██    ██  ██████ ██      ███████     ###
    #>       ###         ██    ██ ████  ████ ██      ██       ██  ██  ██      ██      ██          ###
    #>       ###         ██    ██ ██ ████ ██ █████   ██        ████   ██      ██      █████       ###
    #>       ###         ██    ██ ██  ██  ██ ██      ██         ██    ██      ██      ██          ###
    #>       ###         ██    ██ ██      ██ ███████  ██████    ██     ██████ ███████ ███████     ###
    #>       ########################################################################################
    #> [1] "Starting TimeCycle"
    #> [1] "Pre-Processing Data"
    #> [1] "Computing Periods"
    #> [1] "Pre-Processing Null Distribution"
    #> [1] "Computing Null Distribution"
    #> [1] "Computing Persistence Scores"
    #> [1] "Calculating p-values"
    #> [1] "TimeCycle Completed"
    #> [1] "Analysis Time: 00:00:53"

Check the output and filter for the genes of interest. In this case we
are looking for genes with a period of oscillation between 22 and 26
hours and an FDR  &lt; 0.05.

    library(tidyverse)

    TimeCycleResults %>%
      filter(22 < Period.in.Hours & Period.in.Hours < 26) %>%
      filter(pVals.adj < 0.05) %>%
      summarize(n())
    #>   n()
    #> 1   0

Read `vignette("TimeCycle")` for more details and to learn how to
customise parameters for your time-series data set.


TimeCycle
=========

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![](https://img.shields.io/badge/doi-10.1101/2020.11.19.389981-yellow.svg)](https://doi.org/10.1101/2020.11.19.389981)
[![R build
status](https://github.com/nesscoder/TimeCycle/workflows/R-CMD-check/badge.svg)](https://github.com/nesscoder/TimeCycle/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/badge/devel%20version-1.0.0-blue.svg)](https://github.com/nesscoder/TimeCycle)
[![](https://img.shields.io/github/languages/code-size/nesscoder/TimeCycle.svg)](https://github.com/nesscoder/TimeCycle)

<!-- badges: end -->
<!--
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/TimeCycle.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/TimeCycle)
Add Doi when published

https://doi.org/10.1101/2020.11.19.389981
-->

<img src="man/figures/logo.png" width="200" align="right"/>

> TimeCycle is designed to detect rhythmic genes in circadian
> transcriptomic time-series data. Based on [topological data
> analysis](https://nesscoder.github.io/TimeCycle/articles/TimeCycle.html#takens-theorem-1),
> TimeCycle provides a reliable and efficent reference-free framework
> for cycle detection — handling custom sampling schemes, replicates,
> and missing data.

-   To learn more about the
    [theory](https://nesscoder.github.io/TimeCycle/articles/TimeCycle.html#theory-1)
    and
    [usage](https://nesscoder.github.io/TimeCycle/articles/TimeCycle.html#usage-1)
    of TimeCycle, see our [video
    overview](https://nesscoder.github.io/TimeCycle/articles/TimeCycle.html#video-overview-1)
    and following `vignette("TimeCycle")`.
-   For a comprehensive analysis and discussion of TimeCycle’s
    performance in detecting rhythmic genes, see the accompanying
    [paper](https://doi.org/10.1101/2020.11.19.389981).
-   For details pertaining to the data and source code used in the
    analysis, see the
    [`nesscoder/TimeCycle-data`](https://github.com/nesscoder/TimeCycle-data)
    git repository.

------------------------------------------------------------------------

Installation
------------

TimeCycle has not yet been published on Bioconductor. In the interim
download the development version from GitHub.

<!--
```r
# Install release version from CRAN
install.packages("TimeCycle")

# Install development version from GitHub
devtools::install_github("nesscoder/TimeCycle")
```
-->

    # Install development version from GitHub
    devtools::install_github("nesscoder/TimeCycle")

------------------------------------------------------------------------

Usage
-----

Circadian cycling detection can easily be achieved using TimeCycle’s
main function. Get started with `TimeCycle()` by defining the data and
replicate labels.

In this example, we will use the `zhang2014()` gene expression set
consists of mouse livers sampled every 2-h for 48-h with a single
replicate (i.e. 24 time points). TimeCycle assumes a default circadian
period of 24-h.

-   See `zhang2014()` for additional information about the example data
    set.
-   See [replicate
    labels](https://nesscoder.github.io/TimeCycle/articles/TimeCycle.html#replicate-labels)
    for additional information regarding `repLabel`.

<!-- -->

    library(TimeCycle)

    #set seed for reproducibility with random variables in example usage
    set.seed(1234) 

    TimeCycleResults <- TimeCycle(data = zhang2014, repLabel = rep(1,24))
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
    #> [1] "Analysis Time: 00:00:56"

Once TimeCycle has finished processing, simply check the output and
filter for the genes of interest. In this example, we filter for genes
with a **period of oscillation between 22 and 26** hours and an **FDR
&lt; 0.05**.

    library(tidyverse)

    TimeCycleResults %>%
      filter(22 < Period.in.Hours & Period.in.Hours < 26) %>%
      filter(pVals.adj < 0.05) %>%
      glimpse()
    #> Rows: 1,514
    #> Columns: 7
    #> $ sampleNames     <fct> 1700001C19Rik, 1700010I14Rik, 1700030K09Rik, 1810030O…
    #> $ perScore        <dbl> 0.1183563, 0.1922202, 0.1786031, 0.1501900, 0.1348346…
    #> $ pVals           <dbl> 0.0078, 0.0005, 0.0007, 0.0028, 0.0042, 0.0007, 0.004…
    #> $ pVals.adj       <dbl> 0.04910976, 0.01256641, 0.01274059, 0.02827061, 0.033…
    #> $ Period.in.Hours <dbl> 25.40, 23.50, 24.93, 24.30, 25.73, 24.80, 22.67, 23.3…
    #> $ Amp             <dbl> 0.20, 0.13, 0.07, 0.25, 0.26, 0.16, 0.25, 0.07, 0.13,…
    #> $ Phase.in.Hours  <dbl> 8.5, 7.0, 2.7, 6.5, 2.3, 4.8, 6.2, 4.5, 4.9, 15.8, 5.…

See `vignette("TimeCycle")` for a detailed description of algorithm
design and suggestions for custom parameter selection.

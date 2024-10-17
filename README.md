
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VEMIRT: Variational Expectation Maximization for High-dimensional IRT Models <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The goal of `VEMIRT` is to provide computationally efficient tools for
high-dimensional data, large sample sizes, large item banks, and complex
designs. The package contains several example datasets and functions for

- Parallel analysis
- Exploratory and confirmatory M2PL analysis
- Exploratory and confirmatory M3PL analysis
- Standard error estimates for confirmatory M2PL analysis
- Bootstrap sampling and importance sampling correction for M2PL
  analysis
- Differential item functioning detection for confirmatory M2PL analysis

## Installation

To install this package from source:

1)  **Windows** users may need to install the
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and include
    the checkbox option of installing Rtools to their path for easier
    command line usage. **Mac** users may need to install Xcode command
    line tools by `sudo xcode-select --install` in the terminal, and
    then install [GNU Fortran
    compiler](https://mac.r-project.org/tools/). Most **Linux**
    distributions should already have up to date compilers (or if not
    they can be installed/updated easily).

2)  Install the `devtools` package (if necessary), and install the
    package from [GitHub](https://github.com/) with

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("MAP-LAB-UW/VEMIRT", build_vignettes = T)
torch::install_torch()
```

## Tutorial

After installing the `VEMIRT` package, open its tutorial by running

``` r
vignette("VEMIRT")
```

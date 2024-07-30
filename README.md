
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

1)  Windows users may need to install the
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and include
    the checkbox option of installing Rtools to their path for easier
    command line usage. Mac users will have to download the necessary
    tools from the
    [Xcode](https://apps.apple.com/ca/app/xcode/id497799835?mt=12) and
    its related command line tools (found within Xcode’s Preference Pane
    under Downloads/Components); most Linux distributions should already
    have up to date compilers (or if not they can be updated easily).

2)  Install the `devtools` package (if necessary), and install the
    package from [GitHub](https://github.com/) with

``` r
devtools::install_github("MAP-LAB-UW/VEMIRT")
torch::install_torch()
```

## Tutorial

After installing the `VEMIRT` package, open its tutorial by running

``` r
vignette("VEMIRT")
```

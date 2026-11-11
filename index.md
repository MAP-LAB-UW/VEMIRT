# VEMIRT: Variational Expectation Maximization for High-dimensional IRT Models ![](reference/figures/logo.png)

The goal of `VEMIRT` is to provide computationally efficient tools for
high-dimensional data, large sample sizes, large item banks, and complex
designs. The package contains several example datasets and functions for

- Parallel analysis
- Exploratory and confirmatory M2PL analysis
- Exploratory and confirmatory M3PL analysis
- Standard error estimates for confirmatory M2PL analysis
- Bootstrap sampling and importance sampling correction for M2PL
  analysis
- Differential item functioning detection for confirmatory M1PL and M2PL
  analysis

## Installation

To install this package from source:

1.  **Windows** users may need to install
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) compatible
    with their R version. **macOS** users may need to install the Xcode
    Command Line Tools by running `xcode-select --install` in Terminal,
    followed by the appropriate [GNU Fortran
    compiler](https://mac.r-project.org/tools/). Note that GNU Fortran
    may need to be reinstalled after a macOS upgrade. On **Linux**, the
    required compilers may already be installed; otherwise, they can
    typically be installed or updated through the distribution’s package
    manager.

2.  Install the `pak` package (if necessary), and install the package
    from [GitHub](https://github.com/) with

``` r

if (!require(pak)) install.packages("pak")
options(pkg.build_vignettes = TRUE)
pak::pak("MAP-LAB-UW/VEMIRT", dependencies = TRUE)
torch::install_torch()
```

## Tutorial

After installing the `VEMIRT` package, open its tutorial by running

``` r

vignette("VEMIRT")
```

## Shiny App

After installing the `VEMIRT` package, open its shiny app by running

``` r

VEMIRT::shinyVEMIRT()
```

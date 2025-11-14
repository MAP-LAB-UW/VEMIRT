# Bootstrap Version of GVEM Confirmatory Analysis for M2PL

A bootstrap version of GVEM (i.e., GVEM-BS) can be implemented to
correct the bias on item parameters and compute standard errors under
confirmatory M2PL models

## Usage

``` r
C2PL_bs(gvem_result, boots = 5)
```

## Arguments

- gvem_result:

  a list that includes exploratory or confirmatory GVEM results for M2PL
  models.

- boots:

  the number of bootstrap samples; default is 5

## Value

a list containing the following objects:

- boots_a:

  item discrimination parameters corrected by bootstrap sampling, a \\J
  \times K\\ `matrix`

- boots_b:

  item difficulty parameters corrected by bootstrap sampling, a vector
  of length \\J\\

- sd_a:

  stardard errors of item discrimination parameters, a \\J \times K\\
  matrix

- sd_b:

  stardard errors of item difficulty parameters, a vector of length
  \\J\\

## See also

[`C2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_gvem.md),
[`C2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
gvem_result <- with(C2PL_data, C2PL_gvem(data, model))
C2PL_bs(gvem_result, boots=10)} # }
```

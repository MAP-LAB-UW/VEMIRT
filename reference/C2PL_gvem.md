# Confirmatory M2PL Analysis

Confirmatory M2PL Analysis

## Usage

``` r
C2PL_gvem(u, indic, max.iter = 5000, SE.est = FALSE)
```

## Arguments

- u:

  an \\N \times J\\ `matrix` or a `data.frame` that consists of binary
  responses of \\N\\ individuals to \\J\\ items. The missing values are
  coded as `NA`

- indic:

  a \\J \times K\\ `matrix` or a `data.frame` that describes the factor
  loading structure of \\J\\ items to \\K\\ factors. It consists of
  binary values where 0 refers to the item is irrelevant with this
  factor, 1 otherwise

- max.iter:

  the maximum number of iterations for the EM cycle; default is 5000

- SE.est:

  whether to estimate SE for item parameters using the updated
  supplemented expectation maximization (USEM); default is FALSE

## Value

a list containing the following objects:

- ra:

  item discrimination parameters, a \\J \times K\\ `matrix`

- rb:

  item difficulty parameters, vector of length \\J\\

- reta:

  variational parameters \\\eta(\xi)\\, a \\N \times J\\ matrix

- reps:

  variational parameters \\\xi\\, a \\N \times J\\ matrix

- rsigma:

  population variance-covariance matrix, a \\K \times K\\ matrix

- mu_i:

  mean parameter for each person, a \\K \times N\\ matrix

- sig_i:

  covariance matrix for each person, a \\K \times K \times N\\ array

- n:

  the number of iterations for the EM cycle

- Q_mat:

  factor loading structure, a \\J \times K\\ matrix

- GIC:

  model fit index

- AIC:

  model fit index

- BIC:

  model fit index

- SE:

  Standard errors of item parameters, a \\J \times (K+1)\\ matrix where
  the last column includes SE estimates for item difficulty parameters

## See also

[`C3PL_sgvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C3PL_sgvem.md),
[`C2PL_bs`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_bs.md),
[`C2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(C2PL_data, C2PL_gvem(data, model))} # }
```

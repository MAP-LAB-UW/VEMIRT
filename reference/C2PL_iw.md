# Importance Weighted Version of GVEM Analysis for M2PL Models

An importance weighted version of GVEM (i.e., IW-GVEM) can be
implemented to correct the bias on item parameters under M2PL models

## Usage

``` r
C2PL_iw(u, gvem_result, S = 10, M = 10, max.iter = 10)

E2PL_iw(u, gvem_result, S = 10, M = 10, max.iter = 10)
```

## Arguments

- u:

  a \\N \times J\\ `matrix` or a `data.frame` that consists of binary
  responses of \\N\\ individuals to \\J\\ items. The missing values are
  coded as `NA`

- gvem_result:

  a list that includes exploratory or confirmatory GVEM results for M2PL
  models.

- S:

  the number of times to draw samples;default is 10

- M:

  the number of samples drawn from the variational distributions;default
  is 10

- max.iter:

  the maximum number of iterations for the EM cycle; default is 10

## Value

a list containing the following objects:

- ra:

  item discrimination parameters estimated by GVEM, a \\J \times K\\
  `matrix`

- rb:

  item difficulty parameters estimated by GVEM, vector of length \\J\\

- reta:

  variational parameters \\\eta(\xi)\\, a \\N \times J\\ matrix

- reps:

  variational parameters \\\xi\\, a \\N \times J\\ matrix

- rsigma:

  population variance-covariance matrix estimated by GVEM, a \\K \times
  K\\ matrix

- mu_i:

  mean parameter for each person, a \\K \times N\\ matrix

- sig_i:

  covariance matrix for each person, a \\K \times K \times N\\ array

- n:

  the number of iterations for the EM cycle

- rk:

  factor loadings, a \\J \times K\\ `matrix`, for exploratory analysis
  only

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
  the last column includes SE estimates for item difficulty parameters,
  for confirmatory analysis only

- ur_a:

  item discrimination parameters before conducting the rotation, a \\J
  \times K\\ `matrix`, for exploratory analysis only

- new_a:

  item discrimination parameters estimated by IW-GVEM, a \\J \times K\\
  `matrix`

- new_b:

  item difficulty parameters estimated by IW-GVEM, vector of length
  \\J\\

- new_Sigma_theta:

  population variance-covariance matrix estimated by IW-GVEM, a \\K
  \times K\\ matrix

- best_lr:

  The learning rate used for importance sampling

- best_lb:

  The lower bound value for importance sampling

## See also

[`C2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_gvem.md),
[`E2PL_gvem_rot`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_rot.md),
[`C2PL_bs`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_bs.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
CFA_result <- with(C2PL_data, C2PL_gvem(data, model))
C2PL_iw(C2PL_data$data, CFA_result)} # }
if (FALSE) { # \dontrun{
EFA_result <- with(E2PL_data_C1, E2PL_gvem_lasso(data, model, constrain = constrain, non_pen = non_pen))
E2PL_iw(E2PL_data_C1$data, EFA_result)} # }
```

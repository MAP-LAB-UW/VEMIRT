# Stochastic GVEM for Exploratory M3PL Analysis

Stochastic GVEM for Exploratory M3PL Analysis

## Usage

``` r
E3PL_sgvem_rot(
  u,
  domain,
  samp = 50,
  forgetrate = 0.51,
  mu_b,
  sigma2_b,
  Alpha,
  Beta,
  max.iter = 5000,
  rot = "Promax"
)
```

## Arguments

- u:

  an \\N \times J\\ `matrix` or a `data.frame` that consists of binary
  responses of \\N\\ individuals to \\J\\ items. The missing values are
  coded as `NA`

- domain:

  the number of factors

- samp:

  a subsample for each iteration; default is 50

- forgetrate:

  the forget rate for the stochastic algorithm. The value should be
  within the range from 0.5 to 1. Default is 0.51

- mu_b:

  the mean parameter for the prior distribution of item difficulty
  parameters

- sigma2_b:

  the variance parameter for the prior distribution of item difficulty
  parameters

- Alpha:

  the \\\alpha\\ parameter for the prior distribution of guessing
  parameters

- Beta:

  the \\\beta\\ parameter for the prior distribution of guessing
  parameters

- max.iter:

  the maximum number of iterations for the EM cycle; default is 5000

- rot:

  the post-hoc rotation method: Promax or CF-Quartimax; default is
  `"Promax"`, but may also be `"cfQ"` for conducting the CF-Quartimax
  rotation

## Value

a list containing the following objects:

- ra:

  item discrimination parameters, a \\J \times K\\ `matrix`

- rb:

  item difficulty parameters, vector of length \\J\\

- rc:

  item guessing parameters, vector of length \\J\\

- rs:

  variational parameters \\s\\, a \\N \times J\\ matrix

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

- rk:

  factor loadings, a \\J \times K\\ `matrix`

- GIC:

  model fit index

- AIC:

  model fit index

- BIC:

  model fit index

- ur_a:

  item discrimination parameters before conducting the rotation, a \\J
  \times K\\ `matrix`

## See also

[`E3PL_sgvem_lasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_lasso.md),
[`E3PL_sgvem_adaptlasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_adaptlasso.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
E3PL_sgvem_rot(E3PL_data_C1$data, 3,samp=50,forgetrate=0.51,
mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,rot="Promax")} # }
```

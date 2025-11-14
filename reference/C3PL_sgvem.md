# Stochastic GVEM for Confirmatory M3PL Analysis

Stochastic GVEM for Confirmatory M3PL Analysis

## Usage

``` r
C3PL_sgvem(
  u,
  indic,
  samp = 50,
  forgetrate = 0.51,
  mu_b,
  sigma2_b,
  Alpha,
  Beta,
  max.iter = 5000
)
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

- GIC:

  model fit index

- AIC:

  model fit index

- BIC:

  model fit index

## References

Cho, A. E., Wang, C., Zhang, X., & Xu, G. (2021). Gaussian variational
estimation for multidimensional item response theory. *British Journal
of Mathematical and Statistical Psychology, 74*, 52-85.

Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational
Estimation for Exploratory Item Factor Analysis. *Psychometrika*.
https://doi.org/10.1007/s11336-022-09874-6

## See also

[`C2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_gvem.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(C3PL_data, C3PL_sgvem(data, model, samp=50, forgetrate=0.51, mu_b=0, sigma2_b=4, Alpha=10, Beta=40))} # }
```

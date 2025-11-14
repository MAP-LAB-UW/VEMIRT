# GVEM Algorithm for the Generalized Partial Credit Model

GVEM Algorithm for the Generalized Partial Credit Model

## Usage

``` r
MGPCM_gvem(
  data,
  model = matrix(1, nrow = J, ncol = 4),
  group = rep(1, nrow(data)),
  iter = 2000,
  eps = 1e-05,
  SE = FALSE,
  verbose = TRUE,
  EFA = FALSE
)
```

## Arguments

- data:

  An \\N\times J\\ matrix of item responses where 0 is the minimal
  partial credit score (missing responses should be coded as `NA`)

- model:

  A \\J\times K\\ matrix of loading indicators (K is the Number of
  latent dimension)(all items load on the only dimension by default)

- iter:

  Maximum number of iterations

- eps:

  Termination criterion on numerical accuracy

- SE:

  Whether to calculate the standard errors

- verbose:

  Whether to show the progress

- EFA:

  Whether to rotate the output

## Value

An object of class `vemirt_DIF`, which is a list containing the
following elements:

- ...\$Sigma:

  Group-level covariance matrices

\#'

- ...\$MU:

  Person-level posterior mean vectors

- ...\$a:

  Slopes for group 1

- ...\$b:

  Intercepts for group 1

- ...\$ll:

  Estimated lower bound of log-likelihood

## Author

Yijun Cheng \<chengxb@uw.edu\>

## Examples

``` r
with(MGPCM_gvem, MGPCM_gvem(data, model))
#> Error in eval(substitute(expr), data, enclos = parent.frame()): invalid 'envir' argument of type 'closure'
```

# GVEM Algorithm for the Graded Response Model

GVEM Algorithm for the Graded Response Model

## Usage

``` r
MGRM_gvem(
  data,
  model = matrix(1, ncol(data)),
  method = "GVEM",
  iter = 200,
  tol = 1e-04,
  S = 10,
  M = 10,
  MinDim = 0,
  MaxDim = 0,
  verbose = FALSE,
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

- tol:

  Termination criterion on numerical accuracy

- S:

  Sample size for approximating the expected lower bound (`'IWGVEM'`
  only)

- M:

  Sample size for approximating a tighter lower bound (`'IWGVEM'` only)

- MinDim:

  Minimum num of possible dimensions (`'EFA'` only)

- MaxDim:

  Maximum num of possible dimensions (`'EFA'` only)

- verbose:

  Whether to show the progress

- EFA:

  Whether to run EFA or CFA

- criterion:

  Information criterion for model selection, one of `'GIC'`
  (recommended), `'BIC'`, or `'AIC'`

- c:

  Constant for computing GIC

## Value

An object of class `vemirt_DIF`, which is a list containing the
following elements:

- ...\$SIGMA:

  Person-level posterior covariance matrices

- ...\$MU:

  Person-level posterior mean vectors

- ...\$Sigma:

  Group-level covariance matrices

- ...\$Mu:

  Group-level mean vectors

- ...\$ksi1:

  Variational parameter 1

- ...\$ksi2:

  Variational parameter 2

- ...\$dim:

  Num of dimension between latent variables

- ...\$a:

  Slopes

- ...\$b:

  Intercepts

- ...\$n2vlb:

  Bayesian Information Criterion: `-2*ll+l0*log(N)`

- iter:

  Number(s) of iterations for initialization

## Author

Yijun Cheng \<chengxb@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(MGRM_data, MGRM_gvem(data, method = "IWGVEM", model, EFA = FALSE))} # }
```

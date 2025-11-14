# GVEM Algorithms for DIF Detection in M2PL Models

GVEM Algorithms for DIF Detection in M2PL Models

## Usage

``` r
D2PL_gvem(
  data,
  model = matrix(1, ncol(data)),
  group = rep(1, nrow(data)),
  method = "IWGVEMM",
  Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.1, 1, by = 0.1),
  criterion = "GIC",
  iter = 200,
  eps = 0.001,
  c = 1,
  S = 10,
  M = 10,
  lr = 0.1,
  verbose = TRUE
)
```

## Arguments

- data:

  An \\N\times J\\ binary matrix of item responses (missing responses
  should be coded as `NA`)

- model:

  A \\J\times K\\ binary matrix of loading indicators (all items load on
  the only dimension by default)

- group:

  An \\N\\ dimensional vector of group indicators from `1` to `G` (all
  respondents are in the same group by default)

- method:

  Estimation algorithm, one of `'GVEM'` or `'IWGVEMM'`

- Lambda0:

  A vector of `lambda0` values for \\L_1\\ penalty (`lambda` equals
  `sqrt(N) * lambda0`)

- criterion:

  Information criterion for model selection, one of `'GIC'`
  (recommended), `'BIC'`, or `'AIC'`

- iter:

  Maximum number of iterations

- eps:

  Termination criterion on numerical accuracy

- c:

  Constant for computing GIC

- S:

  Sample size for approximating the expected lower bound (`'IWGVEMM'`
  only)

- M:

  Sample size for approximating a tighter lower bound (`'IWGVEMM'` only)

- lr:

  Learning rate for the Adam optimizer (`'IWGVEMM'` only)

- verbose:

  Whether to show the progress

## Value

An object of class `vemirt_DIF`, which is a list containing the
following elements:

- N:

  Number of respondents

- niter0:

  Number(s) of iterations for initialization

- fit:

  The best (with lowest information criterion) model, which is an
  element of `all`

- best:

  The index of `fit` in `all`

- all:

  A list of models which has the same length as `Lambda0`:

- ...\$lambda0:

  Corresponding element in `Lambda0`

- ...\$lambda:

  `sqrt(N) * lambda0`

- ...\$niter:

  Number(s) of iterations

- ...\$SIGMA:

  Person-level posterior covariance matrices

- ...\$MU:

  Person-level posterior mean vectors

- ...\$Sigma:

  Group-level covariance matrices

- ...\$Mu:

  Group-level mean vectors

- ...\$a:

  Slopes for group 1

- ...\$b:

  Intercepts for group 1

- ...\$gamma:

  D2PL parameters for the slopes

- ...\$beta:

  D2PL parameters for the intercepts

- ...\$RMSE:

  Root mean square error of fitted probability of each item for each
  group

- ...\$ll:

  Estimated lower bound of log-likelihood

- ...\$l0:

  Number of nonzero D2PL parameters in `gamma` and `beta`

- ...\$AIC:

  Akaike Information Criterion: `-2*ll+l0*2`

- ...\$BIC:

  Bayesian Information Criterion: `-2*ll+l0*log(N)`

- ...\$GIC:

  Generalized Information Criterion: `-2*ll+c*l0*log(N)*log(log(N))`

## See also

[`D2PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_pair_em.md),
[`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md),
[`D2PL_lrt`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_lrt.md),
[`coef.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/coef.vemirt_DIF.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(D2PL_data, D2PL_gvem(data, model, group))} # }
```

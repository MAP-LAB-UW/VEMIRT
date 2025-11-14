# EM Algorithms for DIF Detection in M1PL Models

EM Algorithms for DIF Detection in M1PL Models

## Usage

``` r
D1PL_em(
  data,
  model = matrix(1, ncol(data)),
  group = rep(1, nrow(data)),
  a = 1,
  method = "EMM",
  Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.1, 1, by = 0.1),
  level = 10,
  criterion = "BIC",
  iter = 200,
  eps = 0.001,
  c = 1,
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

- a:

  A scalar indicating the common discrimination parameter for all the
  dimensions of all the items (takes `1` by default)

- method:

  Estimation algorithm, one of `'EM'` or `'EMM'`

- Lambda0:

  A vector of `lambda0` values for \\L_1\\ penalty (`lambda` equals
  `sqrt(N) * lambda0`)

- level:

  Accuracy level, either a number for `mvQuad` or a vector indicating
  the grid for each latent dimension

- criterion:

  Information criterion for model selection, one of `'BIC'`
  (recommended), `'AIC'`, or `'GIC'`

- iter:

  Maximum number of iterations

- eps:

  Termination criterion on numerical accuracy

- c:

  Constant for computing GIC

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

- ...\$Sigma:

  Group-level covariance matrices

- ...\$Mu:

  Group-level mean vectors

- ...\$a:

  Slopes for group 1

- ...\$b:

  Intercepts for group 1

- ...\$gamma:

  D1PL parameters for the slopes (all elements are zero)

- ...\$beta:

  D1PL parameters for the intercepts

- ...\$ll:

  Log-likelihood

- ...\$l0:

  Number of nonzero D1PL parameters in `gamma` and `beta`

- ...\$AIC:

  Akaike Information Criterion: `-2*ll+l0*2`

- ...\$BIC:

  Bayesian Information Criterion: `-2*ll+l0*log(N)`

- ...\$GIC:

  Generalized Information Criterion: `-2*ll+c*l0*log(N)*log(log(N))`

## See also

[`D1PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_gvem.md),
[`coef.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/coef.vemirt_DIF.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(D1PL_data, D1PL_em(data, model, group))} # }
```

# EM Algorithm with ADMM for DIF Detection Using Group Pairwise Truncated \\L_1\\ Penalty in 2PL Models

EM Algorithm with ADMM for DIF Detection Using Group Pairwise Truncated
\\L_1\\ Penalty in 2PL Models

## Usage

``` r
D2PL_pair_em(
  data,
  group = rep(1, nrow(data)),
  Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.5, 1.5, by = 0.1),
  Tau = if (length(unique(group)) == 1) 0 else c(Inf, seq(0.05, 0.3, by = 0.05)),
  rho0 = 0.5,
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

- group:

  An \\N\\ dimensional vector of group indicators from `1` to `G` (all
  respondents are in the same group by default)

- Lambda0:

  A vector of `lambda0` values for truncated \\L_1\\ penalty (`lambda`
  equals `sqrt(N) / G * lambda0`)

- Tau:

  A vector of `tau` values for truncated \\L_1\\ penalty (becomes
  \\L_1\\ penalty when `tau` equals `Inf`)

- rho0:

  A value of `rho` for augmented Lagrangian in ADMM (`tau` equals
  `sqrt(N) / G * tau0`)

- level:

  Accuracy level of Gaussian quadrature for `mvQuad`

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

  `sqrt(N) / G * lambda0`

- ...\$tau:

  Corresponding element in `Tau`

- ...\$rho0:

  Same as `rho0` in input

- ...\$rho:

  `sqrt(N) / G * rho0`

- ...\$niter:

  Number(s) of iterations

- ...\$Sigma:

  Group-level covariance matrices

- ...\$Mu:

  Group-level mean vectors

- ...\$a:

  Slopes

- ...\$b:

  Intercepts

- ...\$d.a:

  Group pairwise differences of slopes

- ...\$d.b:

  Group pairwise differences of intercepts

- ...\$u.a:

  Lagrangian multipliers of corresponding elements in `d.a`

- ...\$u.b:

  Lagrangian multipliers of corresponding elements in `d.b`

- ...\$ll:

  Log-likelihood

- ...\$l0:

  Number of nonzero D2PL parameters in `gamma` and `beta`

- ...\$AIC:

  Akaike Information Criterion: `-2*ll+l0*2`

- ...\$BIC:

  Bayesian Information Criterion: `-2*ll+l0*log(N)`

- ...\$GIC:

  Generalized Information Criterion: `-2*ll+c*l0*log(N)*log(log(N))`

## See also

[`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md),
[`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md),
[`D2PL_lrt`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_lrt.md),
[`coef.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/coef.vemirt_DIF.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(D2PL_data, D2PL_pair_em(data, group, Tau = c(Inf, seq(0.01, 0.05, by = 0.01))))} # }
```

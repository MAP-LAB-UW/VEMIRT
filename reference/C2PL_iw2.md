# IW-GVEM Algorithm for Confirmatory M2PL Analysis

IW-GVEM Algorithm for Confirmatory M2PL Analysis

## Usage

``` r
C2PL_iw2(
  data,
  model = matrix(1, ncol(data)),
  criterion = "BIC",
  iter = 200,
  eps = 0.001,
  c = 1,
  S = 10,
  M = 10,
  lr = 0.1,
  SE.level = NULL
)
```

## Arguments

- data:

  An \\N\times J\\ binary matrix of item responses (missing responses
  should be coded as `NA`)

- model:

  A \\J\times K\\ binary matrix of loading indicators (all items load on
  the only dimension by default)

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

  Sample size for approximating the expected lower bound

- M:

  Sample size for approximating a tighter lower bound

- lr:

  Learning rate for the Adam optimizer

- SE.level:

  Accuracy level of Gaussian quadrature for `mvQuad` to compute standard
  errors (SEs are not computed if `SE.level` is `NULL`)

## Value

An object of class `vemirt_DIF`, which is a list containing the
following elements:

- N:

  Number of respondents

- niter0:

  Number(s) of iterations for initialization

- fit:

  The only element of `all`

- best:

  Equal to `1`

- all:

  A list of model which has one element:

- ...\$niter:

  Number(s) of iterations

- ...\$SIGMA:

  Person-level posterior covariance matrices

- ...\$MU:

  Person-level posterior mean vectors

- ...\$Sigma:

  Population covariance matrix

- ...\$Mu:

  Population mean vector

- ...\$a:

  Slopes

- ...\$b:

  Intercepts

- ...\$SE.a:

  Standard errors of `a`

- ...\$SE.b:

  Standard errors of `b`

- ...\$ll:

  Estimated lower bound of log-likelihood

- ...\$l0:

  Number of nonzero elements in `model`

- ...\$AIC:

  Akaike Information Criterion: `-2*ll+l0*2`

- ...\$BIC:

  Bayesian Information Criterion: `-2*ll+l0*log(N)`

- ...\$GIC:

  Generalized Information Criterion: `-2*ll+c*l0*log(N)*log(log(N))`

## See also

[`C2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_gvem.md),
[`C2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md),
[`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md),
[`coef.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/coef.vemirt_DIF.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(C2PL_data, C2PL_iw2(data, model, SE = TRUE))} # }
```

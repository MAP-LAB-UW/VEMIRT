# Stochastic GVEM with Adaptive Lasso Penalty for Exploratory M3PL Analysis

Stochastic GVEM with Adaptive Lasso Penalty for Exploratory M3PL
Analysis

## Usage

``` r
E3PL_sgvem_adaptlasso(
  u,
  indic,
  samp = 50,
  forgetrate = 0.51,
  mu_b,
  sigma2_b,
  Alpha,
  Beta,
  max.iter = 5000,
  constrain = "C1",
  non_pen = NULL,
  gamma = 2
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
  factor, 1 otherwise. For exploratory factor analysis with lasso
  penalty, `indic` should be imposed certain constraints on the a \\K
  \times K\\ sub-matrix to ensure identifiability. The remaining parts
  do not assume any pre-specified zero structure but instead, the
  appropriate lasso penalty would recover the true zero structure. Also
  see `constrain`

- samp:

  a subsample for each iteration; default is 50

- forgetrate:

  the forget rate for the stochastic algorithm. The value should be
  within the range from 0.5 to 1. Default is 0.51

- mu_b:

  the mean parameter for the normal prior distribution of item
  difficulty parameters

- sigma2_b:

  the variance parameter for the normal prior distribution of item
  difficulty parameters

- Alpha:

  the \\\alpha\\ parameter for the beta prior distribution of guessing
  parameters

- Beta:

  the \\\beta\\ parameter for the beta prior distribution of guessing
  parameters

- max.iter:

  the maximum number of iterations for the EM cycle; default is 5000

- constrain:

  the constraint setting: `"C1"` or `"C2"`. To ensure identifiablity,
  `"C1"` sets a \\K \times K\\ sub-matrix of `indic` to be an identity
  matrix.This constraint anchor \\K\\ factors by designating \\K\\ items
  that load solely on each factor respectively. Note that the \\K \times
  K\\ matrix does not have to appear at the top of the `indic` matrix.
  `"C2"` sets the \\K \times K\\ sub-matrix to be a lower triangular
  matrix with the diagonal being ones. That is, there are test items
  associated with each factor for sure and they may be associated with
  other factors as well. Nonzero entries (in the lower triangular part)
  except for the diagonal entries of the sub-matrix are penalized during
  the estimation procedure. For instance, assume \\K=3\\, then the
  `"C2"` constraint will imply the following submatrix:
  \\C2=\begin{bmatrix} 1 & 0 & 0\\ 1 & 1 & 0\\ 1 & 1 &
  1\\\end{bmatrix}\\. As shown, item 1 is allowed to only load on the
  first factor, item 2 will for sure load on the second factor but it
  may also load on the first factor (hence a penalty is added on the
  \\(2,1)\\ element of `"C2"`, i.e., \\C2\_{2,1}\\ ). Item 3 will for
  sure load on the third factor but it may also load on the first two
  factors. However, note that for all remaining items their loading
  vector will all be \\(1, 1, 1)\\ hence indistinguishable from the
  third anchor item. Therefore, we need to alert the algorithm that this
  third anchor item will for sure load on the third factor, and and
  whether or not it loads on the first two factors depends on the
  regularization results. Therefore, we need to specify `"non_pen="` to
  identify the \\K\\th anchor item. Although, `"C2"` is much weaker than
  `"C1"`, it still ensures empirical identifiability. Default is `"C1"`.
  During estimation, under both the `"C1"` and `"C2"` constraints, the
  population means and variances are constrained to be 0 and 1,
  respectively.

- non_pen:

  the index of an item which is associated with each factor to satisfy
  `"C2"`. For `C1`, the input can be `NULL`

- gamma:

  a numerical value of adaptive lasso parameter. Zou (2006) recommended
  three values, 0.5, 1, and 2. The default value is 2.

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

- lbd:

  numerical value of lasso penalty parameter \\\lambda\\

## References

Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational
Estimation for Exploratory Item Factor Analysis. *Psychometrika*.
https://doi.org/10.1007/s11336-022-09874-6

Zou, H. (2006). The adaptive LASSO and its oracle properties. *Journal
of the American Statistical Association, 7*, 1011418–1429.

## See also

[`E3PL_sgvem_rot`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_rot.md),
[`E3PL_sgvem_lasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_lasso.md),
`exampleIndic_efa3pl_c1`, `exampleIndic_efa3pl_c2`

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(E3PL_data_C1, E3PL_sgvem_adaptlasso(data, model,samp=50,forgetrate=0.51,mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,constrain=constrain,non_pen=non_pen,gamma=2))
with(E3PL_data_C2, E3PL_sgvem_adaptlasso(data, model,samp=50,forgetrate=0.51,mu_b=0,sigma2_b=4,Alpha=10,Beta=40,max.iter=5000,constrain=constrain,non_pen=non_pen,gamma=2))} # }
```

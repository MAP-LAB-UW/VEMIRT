# Exploratory M2PL Analysis with Post-hoc Rotation

Exploratory M2PL Analysis with Post-hoc Rotation

## Usage

``` r
E2PL_gvem_rot(u, domain, max.iter = 5000, rot = "Promax")
```

## Arguments

- u:

  an \\N \times J\\ `matrix` or a `data.frame` that consists of binary
  responses of \\N\\ individuals to \\J\\ items. The missing values are
  coded as `NA`

- domain:

  the number of factors

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

- rk:

  factor loadings, a \\J \times K\\ `matrix`

- Q_mat:

  factor loading structure, a \\J \times K\\ matrix

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

[`E2PL_gvem_lasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_lasso.md),
[`E2PL_gvem_adaptlasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_adaptlasso.md)

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
E2PL_gvem_rot(E2PL_data_C1$data, domain=5,max.iter=3000)
E2PL_gvem_rot(E2PL_data_C1$data, domain=5,rot="cfQ")} # }
```

# Likelihood Ratio Test for DIF Detection in M2PL Models

Likelihood Ratio Test for DIF Detection in M2PL Models

## Usage

``` r
D2PL_lrt(data, model, group, unif = F)
```

## Arguments

- data:

  An \\N\times J\\ binary matrix of item responses

- model:

  A \\J\times K\\ binary matrix of loading indicators

- group:

  An \\N\\ dimensional vector of group indicators (integers from `1` to
  `G`)

- unif:

  Whether to detect uniform D2PL only

## Value

A list:

- Sigma:

  Group-level posterior covariance matrices

- Mu:

  Group-level posterior mean vectors

- a:

  Slopes for group 1

- b:

  Intercepts for group 1

- gamma:

  D2PL parameters for the slopes

- beta:

  D2PL parameters for the intercepts

## See also

[`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md),
[`D2PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_pair_em.md),
[`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md)

## Author

Ruoyi Zhu \<zhux0445@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
with(D2PL_data, D2PL_lrt(data, model, group))} # }
```

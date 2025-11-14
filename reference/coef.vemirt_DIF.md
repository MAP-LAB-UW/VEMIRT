# Extract Parameter Estimates from DIF 2PL Analysis

Extract Parameter Estimates from DIF 2PL Analysis

## Usage

``` r
coef(object, criterion = NULL)
```

## Arguments

- object:

  An object of class `vemirt_DIF`

- criterion:

  Information criterion for model selection, one of `'AIC'`, `'BIC'`,
  `'GIC'`, or the constant for computing GIC, otherwise use the
  criterion specified when fitting the model(s)

## See also

[`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md),
[`D2PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_pair_em.md),
[`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<wlyu4@uw.edu\>

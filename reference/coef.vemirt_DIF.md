# Extract Parameter Estimates from DIF Analysis

Extract Parameter Estimates from DIF Analysis

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

[`D1PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_em.md),
[`D1PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_pair_em.md),
[`D1PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_gvem.md),
[`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md),
[`D2PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_pair_em.md),
[`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md),
[`print.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/print.vemirt_DIF.md),
[`summary.vemirt_DIF`](https://MAP-LAB-UW.github.io/VEMIRT/reference/summary.vemirt_DIF.md)

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

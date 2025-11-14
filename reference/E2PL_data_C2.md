# Simulated Data Set for Exploratory M2PL Analysis Under C2 Constraint

Responses are simulated based on an M2PL model with 3 factors. The true
factor correlations are set as 0.5.

## Usage

``` r
E2PL_data_C2
```

## Format

A list of components of the data set:

|             |                                                          |
|-------------|----------------------------------------------------------|
| ` data`     | Item responses                                           |
|             |                                                          |
| `model`     | Loading indicators for (adaptive) lasso penalty          |
|             |                                                          |
| `constrain` | Constraint for model identification (`'C2'`)             |
|             |                                                          |
| `non_pen`   | Index of an item that is associated with all the factors |
|             |                                                          |
| `params`    | True parameters used for generating the item responses   |

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

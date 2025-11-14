# Simulated Data Set for Exploratory M3PL Analysis Under C1 Constraint

Responses are simulated based on an M3PL model with 3 factors. The true
factor correlations are set as 0.5.

## Usage

``` r
E3PL_data_C1
```

## Format

A list of components of the data set:

|             |                                                                              |
|-------------|------------------------------------------------------------------------------|
| ` data`     | Item responses                                                               |
|             |                                                                              |
| `model`     | Loading indicators for (adaptive) lasso penalty                              |
|             |                                                                              |
| `constrain` | Constraint for model identification (`'C1'`)                                 |
|             |                                                                              |
| `non_pen`   | Index of an item that is associated with all the factors (`NULL` under `C1`) |
|             |                                                                              |
| `params`    | True parameters used for generating the item responses                       |

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

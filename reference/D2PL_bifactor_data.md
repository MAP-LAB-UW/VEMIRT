# Simulated Data Set for DIF M2PL Bifactor Analysis

Simulated Data Set for DIF M2PL Bifactor Analysis

## Usage

``` r
D2PL_bifactor_data
```

## Format

A list of components of the data set:

|              |                                                                                        |
|--------------|----------------------------------------------------------------------------------------|
| ` data`      | Item responses                                                                         |
|              |                                                                                        |
| `model`      | Loading indicators                                                                     |
|              |                                                                                        |
| `ndim`       | The first `ndim` dimensions are substantive factors, while others are nuisance factors |
|              |                                                                                        |
| `impact`     | Whether nuisance factors have impact                                                   |
|              |                                                                                        |
| `group`      | Group indicators                                                                       |
|              |                                                                                        |
| `j`          | Indices of DIF items                                                                   |
|              |                                                                                        |
| `params`     | A list of true parameters used for generating the item responses:                      |
|              |                                                                                        |
| ` ...$a`     | Slopes                                                                                 |
|              |                                                                                        |
| ` ...$b`     | Negated intercepts                                                                     |
|              |                                                                                        |
| ` ...$theta` | Latent traits                                                                          |

## Author

Yijun Cheng \<chengxb@uw.edu\>, Weicong Lyu \<weiconglyu@um.edu.mo\>

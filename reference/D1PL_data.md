# Simulated Data Set for DIF M1PL Analysis

Simulated Data Set for DIF M1PL Analysis

## Usage

``` r
D1PL_data
```

## Format

A list of components of the data set:

|              |                                                                   |
|--------------|-------------------------------------------------------------------|
| ` data`      | Item responses                                                    |
|              |                                                                   |
| `model`      | Loading indicators                                                |
|              |                                                                   |
| `group`      | Group indicators                                                  |
|              |                                                                   |
| `j`          | Number of DIF items (the first `j` items have DIF)                |
|              |                                                                   |
| `params`     | A list of true parameters used for generating the item responses: |
|              |                                                                   |
| ` ...$a`     | Slopes                                                            |
|              |                                                                   |
| ` ...$b`     | Negated intercepts                                                |
|              |                                                                   |
| ` ...$theta` | Latent traits                                                     |

## Author

Weicong Lyu \<weiconglyu@um.edu.mo\>

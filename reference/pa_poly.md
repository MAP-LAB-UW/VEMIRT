# Parallel analysis using polychoric correlation

Identify the number of factors

## Usage

``` r
pa_poly(data, n.iter = 10, figure = TRUE)
```

## Arguments

- data:

  a \\N \times J\\ `matrix` or a `data.frame` that consists of the
  responses of \\N\\ individuals to \\J\\ items without any missing
  values. The responses are binary or polytomous.

- n.iter:

  Number of simulated analyses to perform

- figure:

  By default, `pa_poly` draws an eigenvalue plot. If FALSE, it
  suppresses the graphic output

## Value

`pa_poly` returns a `data.frame` with the eigenvalues for the real data
and the simulated data.

## Author

Jiaying Xiao \<jxiao6@uw.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
pa_poly(C2PL_data$data, n.iter=20)} # }
```

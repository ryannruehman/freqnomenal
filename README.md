
# fourieR

<img src="man/figures/fourieR.png" width="120" />
 
<!-- badges: start -->
<!-- badges: end -->

The goal of fourieR is to compute the Discrete Fourier Transform and Power Spectral Density
    of a signal using real-valued least squares decomposition, avoiding the 
    use of complex numbers.

## Installation

You can install the development version of fourieR from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ryannruehman/fourieR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fourieR)
signal <- c(2, 4, 1, -7, 2, 3, 5, -4)
time <- c(0, 1, 2, 3, 4, 5, 6, 7)
dft(signal, time)
psd(signal, time)
```


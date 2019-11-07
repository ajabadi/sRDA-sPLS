RDA & PLS
================

``` r
suppressMessages({
    library(microbenchmark)
    library(Rfast)
})
```

``` r
sourceCpp("../src/test.cpp")
```

# matrices

``` r
A <- matrix(rnorm(1000000), nrow = 300, ncol = 1000)
```

    ## Warning in matrix(rnorm(1e+06), nrow = 300, ncol = 1000): data length
    ## [1000000] is not a sub-multiple or multiple of the number of rows [300]

``` r
B <- matrix(rnorm(10000), nrow = 300, ncol = 500)
dim(A); dim(B)
```

    ## [1]  300 1000

    ## [1] 300 500

## consistency

``` r
res<-list(
`%*%`=t(A) %*% B,
crossprd=crossprod(A, B),
matmult=mat.mult(t(A), B),
eMM=eigenMatMult(t(A), B),
eMMM=eigenMapMatMult(t(A), B))

## all output the same matrix?
all.identical <- function(x) length(unique(x)) == 1
all.identical(res)
```

    ## [1] TRUE

# run times

``` r
microbenchmark(list = list(
    `%*%`=t(A) %*% B,
    crossprd=crossprod(A, B),
    matmult=mat.mult(t(A), B),
    eMM=eigenMatMult(t(A), B),
    eMMM=eigenMapMatMult(t(A), B)) , times = 10)
```

    ## Unit: nanoseconds
    ##      expr min lq mean median uq max neval
    ##       %*%   4  5 28.8      6  6 238    10
    ##  crossprd   4  5 17.5      5  6 129    10
    ##   matmult   3  4 12.1      5  6  77    10
    ##       eMM   4  4 12.6      5  6  81    10
    ##      eMMM   4  5 22.7      5  6 181    10

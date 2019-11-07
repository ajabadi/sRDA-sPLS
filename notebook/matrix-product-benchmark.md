Efficient Matrix Multiplication in R
================

``` r
suppressMessages({
    library(microbenchmark)
    library(Rfast)
})
```

``` r
sourceCpp("../src/utils/test.cpp")
```

# matrices - small

``` r
colsA <- 1000
colsB <- 2000
A <- matrix(rnorm(300*colsA), nrow = 300, ncol = colsA)
B <- matrix(rnorm(300*colsB), nrow = 300, ncol = colsB)
dim(A); dim(B)
```

    ## [1]  300 1000

    ## [1]  300 2000

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
#> TRUE
```

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
    ##       %*%   5  6 22.8    7.0  7 169    10
    ##  crossprd   5  6 18.9    7.0  7 131    10
    ##   matmult   6  6 46.8    7.0  7 405    10
    ##       eMM   6  6 30.4    7.5  8 240    10
    ##      eMMM   6  6 22.0    6.5  7 161    10

# matrices - larger

``` r
colsA <- 3000
colsB <- 5000
A <- matrix(rnorm(300*colsA), nrow = 300, ncol = colsA)
B <- matrix(rnorm(300*colsB), nrow = 300, ncol = colsB)
dim(A); dim(B)
```

    ## [1]  300 3000

    ## [1]  300 5000

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
#> TRUE
```

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
    ##       %*%   5  6 26.9    7.5 12 192    10
    ##  crossprd   5  7 26.2    9.5 11 181    10
    ##   matmult   6  7 48.9    9.5 12 406    10
    ##       eMM   7 10 28.0   11.0 12 186    10
    ##      eMMM   6  8 77.0    9.5 11 686    10

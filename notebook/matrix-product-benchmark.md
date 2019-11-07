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
A <- matrix(rnorm(300*colsA), nrow = 300, ncol = colsA) %>% t()
B <- matrix(rnorm(300*colsB), nrow = 300, ncol = colsB)
dim(A); dim(B)
```

    ## [1] 1000  300

    ## [1]  300 2000

## consistency

``` r
res<-list(
`%*%`=A %*% B,
# crossprd=crossprod(A, B),
matmult=mat.mult(A, B),
eMM=eigenMatMult(A, B),
eMMM=eigenMapMatMult(A, B))

## all output the same matrix?
all.identical <- function(x) length(unique(x)) == 1
all.identical(res)
#> TRUE
```

# run times

``` r
microbenchmark(
A %*% B,
mat.mult(A, B),
eigenMatMult(A, B),
eigenMapMatMult(A, B) , times = 10)
```

    ## Unit: milliseconds
    ##                   expr       min        lq     mean   median       uq
    ##                A %*% B 370.79607 371.73450 394.3669 375.7986 439.1958
    ##         mat.mult(A, B) 203.36440 210.68277 223.1740 220.5062 240.9514
    ##     eigenMatMult(A, B)  97.21408  99.99507 105.0663 101.0127 111.9635
    ##  eigenMapMatMult(A, B)  96.25533  98.25440 105.8678 102.6354 112.9735
    ##       max neval
    ##  443.8928    10
    ##  243.9406    10
    ##  114.4928    10
    ##  120.5286    10

# matrices - larger

``` r
colsA <- 3000
colsB <- 5000
A <- matrix(rnorm(300*colsA), nrow = 300, ncol = colsA)%>% t()
B <- matrix(rnorm(300*colsB), nrow = 300, ncol = colsB)
dim(A); dim(B)
```

    ## [1] 3000  300

    ## [1]  300 5000

## consistency

``` r
res<-list(
`%*%`=A %*% B,
# crossprd=crossprod(A, B),
matmult=mat.mult(A, B),
eMM=eigenMatMult(A, B),
eMMM=eigenMapMatMult(A, B))

## all output the same matrix?
all.identical <- function(x) length(unique(x)) == 1
all.identical(res)
#> TRUE
```

# run times

``` r
microbenchmark(
A %*% B,
mat.mult(A, B),
eigenMatMult(A, B),
eigenMapMatMult(A, B) , times = 10)
```

    ## Unit: milliseconds
    ##                   expr       min        lq      mean    median        uq
    ##                A %*% B 2947.5140 2986.0156 3018.2000 3017.2993 3055.3269
    ##         mat.mult(A, B) 1624.0355 1628.8686 1663.5077 1645.1228 1698.0027
    ##     eigenMatMult(A, B)  737.6062  743.6640  779.7215  759.7714  776.1901
    ##  eigenMapMatMult(A, B)  735.6342  739.2713  760.9625  749.6793  771.7794
    ##        max neval
    ##  3090.2043    10
    ##  1776.0234    10
    ##   908.4489    10
    ##   843.8530    10

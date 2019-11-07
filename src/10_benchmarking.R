library(microbenchmark)
library(Rcpp)
library(Rfast)

A <- matrix(rnorm(1000000), nrow = 300, ncol = 1000)
B <- matrix(rnorm(10000), nrow = 300, ncol = 500)
dim(A); dim(B)

# sourceCpp("src/test.cpp")

res<-list(
`%*%`=t(A) %*% B,
crossprd=crossprod(A, B),
matmult=mat.mult(t(A), B),
eMM=eigenMatMult(t(A), B),
eMMM=eigenMapMatMult(t(A), B))

all.identical <- function(x) length(unique(x)) == 1

all.identical(res)

microbenchmark(list = list(
    `%*%`=t(A) %*% B,
    crossprd=crossprod(A, B),
    matmult=mat.mult(t(A), B),
    eMM=eigenMatMult(t(A), B),
    eMMM=eigenMapMatMult(t(A), B)) , times = 10)

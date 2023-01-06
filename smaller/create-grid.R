# `SMALL` VERSION
library(Rcpp)
library(RcppArmadillo)
sourceCpp('funs.cpp')
# Creating grid of lambda values, may add log(Z) in future
.f <- function(x) format(round(x, 2), nsmall = 2)
generate.grid <- function(mus, nus){
    M <- structure(poly_mat(mus, nus),
                   dimnames = list(.f(mus), .f(nus)))
    M
}

# Short test
m <- runif(100, 0.98, 1.99);  n <- runif(100, 1.22, 1.65)
M <- generate.grid(seq(min(round(m, 2)), max(round(m, 2)), .01),
                   seq(min(round(n, 2)), max(round(n, 2)), .01))
rm(m,n,M)


# Lookup ------------------------------------------------------------------
get.indices <- function(this, f = rownames){
    this <- round(this, 2);
    parent <- as.numeric(f(M)) # maybe remove
    match(this, parent)
}
grid.lookup <- function(mus, nus){
    # Work out the indices of mus, nus to lookup
    mu.lookup <- get.indices(mus); nu.lookup <- get.indices(nus, colnames)
    M[mu.lookup, nu.lookup]
}

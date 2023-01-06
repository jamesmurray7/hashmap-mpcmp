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

# Lookup ------------------------------------------------------------------
.r <- function(x) round(x, 2)

get.indices <- function(this, f = rownames){
    parent <- as.numeric(f(M)) # maybe remove
    match(this, parent)
}

grid.lookup <- function(mus, nus){
    mus <- unique(.r(mus)); nus <- unique(.r(nus))
    # Work out the indices of mus, nus to lookup
    mu.lookup <- unique(get.indices(mus))
    nu.lookup <- unique(get.indices(nus, colnames))

    out <- M[mu.lookup, nu.lookup]

    # Above matrix generates some <NA> values...
    missing.nu <- which(is.na(nu.lookup)); has.nu <- which(!is.na(nu.lookup))
    missing.mu <- which(is.na(mu.lookup)); has.mu <- which(!is.na(mu.lookup))
    if(length(missing.nu)){
        missing.vals <- poly_mat(mus, nus[missing.nu])
        out[, missing.nu] <- missing.vals
        colnames(out)[missing.nu] <- .f(nus[missing.nu])
    }
    if(length(missing.mu)){
        missing.vals <- poly_mat(mus[missing.mu], nus[has.nu])
        # Need to be careful, nu has already poopulated this mu value.
        out[missing.mu, has.nu] <- missing.vals
        rownames(out)[missing.mu] <- .f(mus[missing.mu])
    }
    out
}

m2 <- runif(100, 0.98, 1.99); n2 <- runif(100, 1.22, 1.65)
M2 <- grid.lookup(m2, n2)
rm(m,m2,n,n2,M)

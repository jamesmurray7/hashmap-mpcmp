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
# m <- runif(100, 0.98, 1.99);  n <- runif(100, 1.22, 1.65)
# M <- generate.grid(seq(min(round(m, 2)), max(round(m, 2)), .01),
#                    seq(min(round(n, 2)), max(round(n, 2)), .01))

# Lookup ------------------------------------------------------------------
.r <- function(x) round(x, 2)

get.indices <- function(this, f = rownames){
    parent <- as.numeric(f(M)) # maybe remove
    match(this, parent)
}

grid.lookup <- function(mus, nus){
    rmus <- .r(mus); rnus <- .r(nus)
    # Work out the indices of mus, nus to lookup
    mu.lookups <- get.indices(rmus)
    nu.lookups <- get.indices(rnus, colnames)

    rtn <- mapply(function(m, n){
        if(!is.na(m * n)){
            return(M[m, n])
        }else{
            return(NA)
        }
    }, m = mu.lookups, n = nu.lookups, SIMPLIFY = T)

    # Populate non-matching ones
    w <- which(is.na(rtn))
    if(length(w)){ # Lazy --> Assumes nu is scalar value (i.e. one for all...)
        rtn[w] <- mapply(function(i) poly_solve(mus[i], nus), i = w, SIMPLIFY = T)

        # Lines to update...
        # R <- as.numeric()
    }
    rtn

    # Some way to add these to grid?
}

# m2 <- runif(100, 0.98, 1.99); n2 <- runif(100, 1.22, 1.65)
# grid.lookup(m2,n2)
# rm(m,m2,n,n2,M)

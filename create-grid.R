library(Rcpp)
library(RcppArmadillo)
sourceCpp('funs.cpp')
# Creating grid of lambda values, may add log(Z) in future
mus <- nus <- seq(0.01, 9.99, 0.01)
M <- poly_mat(mus,nus)
# Set dimnames
# Way one --> '0.45,1.22' == 'mu,nu'.
.f <- function(x) format(round(x, 2), nsmall = 2)
M <- structure(M,
               dimnames = list(.f(mus), .f(nus)))


# Lookup items from the grid ----------------------------------------------
grid.lookup <- function(mu, nu){
    l <- length(mu); k <- length(nu)
    lam <- numeric(l)
    mus <- round(mu, 2)*100; nus <- round(nu, 2) * 100
    for(i in 1:l){
        lam[i] <- M[mus[i], nus[i]]
    }
    lam
}

grid.lookup2 <- function(mu, nu){
    mu <- round(mu, 2)*100; nu <- round(nu, 2) * 100
    M[mu,nu]
}

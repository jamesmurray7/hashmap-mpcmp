library(Rcpp); library(RcppArmadillo)
sourceCpp('funs.cpp')
# Trying out with hashmap -------------------------------------------------
mus <- nus <- seq(0.01, 9.99, 0.01)
M <- poly_mat(mus,nus)
# Set dimnames
# Way one --> '0.45,1.22' == 'mu,nu'.
.f <- function(x) format(round(x, 2), nsmall = 2)
M <- structure(M,
               dimnames = list(.f(mus), .f(nus)))


# Function to create a hashmap --------------------------------------------
library(hashmap)
mat2hash <- function(mat){
    R <- rownames(mat); C <- colnames(mat)
    lookup <- c(mat) # column-by-column concatenation
    keys <- apply(expand.grid(R,C), 1, paste, collapse=',')
    hashmap(keys, lookup)
}
HH <- mat2hash(M)

# Second way --> Looking up a double faster than converting to string first?
mat2hash_b <- function(mus, nus, M){
    keys <- outer(mus, nus, '*')
    hashmap(c(keys), c(M))
}

# Seems to work quite well!
HH[['1.52,1.20']] # poly_solve(1.52, 1.20)
HH[['0.01,0.01']] # poly_solve(0.01, 0.01)
HH[[c("1.00,1.01", "1.92,0.20", "0.20,1.24")]]
poly_solve(1,1.01); poly_solve(1.92,.2); poly_solve(.2,1.24)


# Functions to lookup from grid and hash ----------------------------------
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
    l <- length(mu); k <- length(nu)
    lam <- numeric(l)
    mus <- round(mu, 2)*100; nus <- round(nu, 2) * 100
    M[mus,nus]
}

hash.lookup <- function(mu, nu){
    lookups <- paste0(.f(mu), ',', .f(nu))
    HH[[lookups]]
}
grid.lookup(seq(.55,1,.01), seq(1, 1.45, .01))
hash.lookup(seq(.55,1,.01), seq(1, 1.45, .01)) # arbitrary e.g.


# Benchmarking ------------------------------------------------------------
# Create a scenario where we have N lookups to do (which are random numbers).
library(microbenchmark)
# Function
benchN <- function(N, times = 1000){
    mus <- runif(N, 0.02, 9.99)
    nus <- runif(N, 0.02, 9.99)
    bench <- microbenchmark(
        `grid` = {
            G <- grid.lookup(mus, nus)
        },
        `hash` = {
            H <- hash.lookup(mus, nus)
        },
        times = times
    )
    bench
}
b <- benchN(1000, 1000)


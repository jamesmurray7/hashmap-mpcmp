rm(list=ls())
library(glmmTMB)
library(hashmap)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)

# Loose idea of how may fit using a hashmap.
source('simData.R')
data <- simData1()$data
data$dum <- 1

#### Function to extract some data summaries needed in updates
summarydata <- function(y)
{
    N <- length(y)
    unique <- 0:max(y) #sort(unique(data[, col]))
    lfac <- lfactorial(unique)
    lfac_y <- lfactorial(y)
    sum_lfac_y <- sum(lfactorial(y))
    sum_y <- sum(y)
    list(N = N, unique = unique, lfac = lfac, lfac_y = lfac_y, sum_lfac_y = sum_lfac_y,
         sum_y = sum_y)
}

source('create-grid.R')
source('create-hashmap.R')

fit <- function(simdata, N = 1e3, lookup = 'hashmap'){

    if(lookup == 'hashmap') .lookup <- hash.lookup else .lookup <- grid.lookup2

    f <- suppressWarnings(glmmTMB::glmmTMB(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10 - 1, data = simdata,
                                           disp = ~ 1, family = glmmTMB::compois()))

    X <- getME(f, 'X')
    y <- simdata$y

    # Point estimates
    nu.hat <- exp(-fixef(f)$disp)
    beta.hat <- fixef(f)$cond

    # Variance-covariance of MLEs
    V <- vcov(f,T)
    Sigma.beta <- V[1:10,1:10]

    # Covariance matrix for proposal distribution
    Sigma.beta <- 0.7*Sigma.beta

    # Prior mean and variance for Beta
    p <- dim(X)[2]
    mu0 <- rep(0,p)
    Sigma0 <- diag(p)*10^5

    # Initial value for Beta
    beta0 <- beta.hat
    mu0 <- exp(X %*% beta0)
    if (max(mu0) >= 10){mu0 <- pmin(9.99, mu0)}
    if (min(mu0) < 0.01){mu0 <- pmax(0.01, mu0)}
    # Initial value for nu
    nu0 <- nu.hat

    # Setup
    betas <- matrix(0, nrow = N, ncol = p)
    nus <- matrix(0, nrow=N, ncol = 1)
    colnames(betas) <-  paste0('x',attr(X, 'assign'))
    colnames(nus) <- "dispersion"
    iter_times <- rep(0, N)
    ptm <- proc.time()

    # summary of data
    summ <- summarydata(simdata$y)
    lfac <- summ$lfac
    lfac_y <- summ$lfac_y
    n <- summ$N
    unique <- summ$unique
    sum_y <- summ$sum_y

    for(i in 1:N){
        s <- proc.time()[3]
        if(i %% 100 == 0) cat(sprintf("Iteration %d done\r", i))

        # Beta update
        beta1 <- rmvnorm(1, beta0, Sigma.beta)
        mu1 <- exp(tcrossprod(X, beta1))
        if (max(mu1) >= 10){mu1 <- pmin(9.9, mu1)}
        if (min(mu1) < 0.01){mu1 <- pmax(0.01, mu1)}
        # Grid part

        log_lambda <- log(.lookup(mu0, nu0))
        log_lambda1 <- log(.lookup(mu1, nu0))

        # Calculate ln(G) directly (depends on lambda)
        lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(rep(nu0, n), lfac))))
        lG1 <- log(rowSums(exp(tcrossprod(log_lambda1, unique) - tcrossprod(rep(nu0, n), lfac))))
        denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG))
        laccept <- min(0, denrat)
        accept <- (log(runif(1)) < laccept)
        if (accept){
            beta0 <- beta1
            mu0 <- mu1
        }
        betas[i,] = beta0

        # Dispersion update
        nu1 <- rexp(1, 1/nu0)
        if (nu1 < 0.01){nu1 <- 0.01}
        if (nu1 > 9.9){nu1 <- 9.9}
        log_lambda <- log(.lookup(mu0, nu0))
        log_lambda1 <- log(.lookup(mu0, nu1))
        lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(rep(nu0, n), lfac))))
        lG1 <- log(rowSums(exp(tcrossprod(log_lambda1, unique) - tcrossprod(rep(nu1, n), lfac))))
        denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG) + (nu0 - nu1)*lfac_y)
        nu.ratio <- nu1/nu0
        denrat <- denrat + log(nu.ratio) + nu.ratio - 1/nu.ratio
        laccept <- min(0, denrat)
        accept <- (log(runif(1)) < laccept)
        if (accept){
            nu0 <- nu1
        }
        nus[i,] <- nu0
        iter_times[i] <- proc.time()[3] - s
    }

    paras <- cbind(betas, nus)
    cpu_time <- proc.time() - ptm
    list(cpu_time = cpu_time[3], paras = paras, iter_times = iter_times)
}

# Below implies that as more grid lookups are done, hashmap becomes orders slower.
a <- fit(data)
b <- fit(data, lookup='grid')
cat(sprintf("Hashmap: %.2fs, grid: %.2fs.\n", a$cpu_time, b$cpu_time))

a <- fit(data, N = 5e3)
b <- fit(data, N = 5e3, lookup = 'grid')
cat(sprintf("Hashmap: %.2fs, grid: %.2fs.\n", a$cpu_time, b$cpu_time))
# Is this because the hashmap is so large?
# Carry out some investigation into size of tables --> Does hashmap perform better if it's smaller?

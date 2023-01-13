# `SMALL` version
rm(list=ls())
library(glmmTMB)
library(hashmap)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(mpcmp)

# Loose idea of how may fit using a hashmap.
source('smaller/create-grid.R')
source('smaller/create-hashmap.R')
source('simData.R')

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

simdata <- simData1()$data

fit <- function(simdata, N = 1e3, lookup = 'hashmap', CI = .95, verbose = TRUE, update.grid = T){

    if(lookup == 'hashmap') .lookup <- hash.lookup else .lookup <- function(m,n,u) grid.lookup(m,n,update.grid)

    f <- suppressWarnings(glmmTMB::glmmTMB(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10 - 1, data = simdata,
                                           disp = ~ 1, family = glmmTMB::compois()))

    # f2 <- glm.cmp(y ~ . -1, data = simdata)

    X <- getME(f, 'X')
    y <- simdata$y

    # Point estimates
    nu.hat <- exp(-fixef(f)$disp)
    beta.hat <- fixef(f)$cond

    # Variance-covariance of MLEs
    V <- vcov(f,T)
    Sigma.beta <- V[1:10,1:10]
    Vnu <- V[11,11]

    # Setup grid (and hash) based on smaller subset.
    p <- dim(X)[2]
    upper.beta <- sapply(1:p, function(x) qnorm(1-((1 - CI)/2), beta.hat[x], sqrt(Sigma.beta[x,x])))
    lower.beta <- sapply(1:p, function(x) qnorm((1 - CI)/2, beta.hat[x], sqrt(Sigma.beta[x,x])))
    upper.nu <- qnorm(1-((1 - CI)/2), fixef(f)$disp, sqrt(Vnu))
    lower.nu <- qnorm((1 - CI)/2, fixef(f)$disp, sqrt(Vnu))

    # Work out grid values
    mu.lower <- exp(X %*% lower.beta); mu.upper <- exp(X %*% upper.beta)
    nu.lower <- exp(-upper.nu); nu.upper <- exp(-lower.nu)

    mus.for.grid <- seq(.r(quantile(mu.lower)[2]), .r(quantile(mu.upper)[3]), by = 0.01)
    if (max(mus.for.grid) >= 10){mus.for.grid <- pmin(9.99, mus.for.grid)}
    if (min(mus.for.grid) < 0.01){mus.for.grid <- pmax(0.01, mus.for.grid)}
    nus.for.grid <- seq(.r(nu.lower), .r(nu.upper), by = 0.01)
    if(verbose){
        cat(sprintf("Based on the %d%% confidencce interval for values of beta and nu,\n", 100 * CI))
    #    cat(sprintf("Creating a %d x %d grid for mu x nu values.\n", length(mus.for.grid), length(nus.for.grid)))
    }

    M <<- generate.grid(mus.for.grid, nus.for.grid)
    if(lookup == 'hashmap'){
        HH <<- mat2hash(M)
        if(verbose)
            cat(sprintf("And created hashmap of size %d.\n", HH$size()))
    }

    # Covariance matrix for proposal distribution
    Sigma.beta <- 0.7*Sigma.beta

    # Prior mean and variance for Beta
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
        if(lookup == 'hashmap'){
            if(i %% 100 == 0) cat(sprintf("Iteration %d done, Hashmap size is now %d\r", i, HH$size()))
        }else{
            if(i %% 100 == 0) cat(sprintf("Iteration %d done, grid is now %d x %d\r", i, nrow(M), ncol(M)))
        }
            

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
a <- fit(simdata, CI = .3)
b <- fit(simdata, lookup='grid', CI = .3, update.grid = T)
cat(sprintf("Hashmap: %.2fs, grid: %.2fs.\n", a$cpu_time, b$cpu_time))

# A lot more draws
# NB hash seems faster, this b/c grid doesn _not_ update if {nu, nu} not
# found (I think).
a2 <- fit(data, CI = .5, N = 5e3)
b2 <- fit(data, CI = .5, N = 5e3, lookup = 'grid')

cat(sprintf("Hashmap: %.2fs, grid: %.2fs.\n", a2$cpu_time, b2$cpu_time))
apply(a2$paras, 2, mean); apply(a2$paras, 2, sd)
apply(b2$paras, 2, mean); apply(b2$paras, 2, sd)

# Create more data-points
data2 <- simData1(n = 1000)
data2 <- data2$data
a2b <- fit(data2, CI = .5, N = 5e3)
b2b <- fit(data2, CI = .5, N = 5e3, lookup = 'grid')
cat(sprintf("Hashmap: %.2fs, grid: %.2fs.\n", a2b$cpu_time, b2b$cpu_time))

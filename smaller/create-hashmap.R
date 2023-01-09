# `SMALL` version
library(hashmap)
.f <- function(x) format(round(x, 2), nsmall = 2)
# Creating hashmap --------------------------------------------------------
mat2hash <- function(mat){
    R <- rownames(mat); C <- colnames(mat)
    lookup <- c(mat) # column-by-column concatenation
    keys <- apply(expand.grid(R,C), 1, paste, collapse=',')
    hashmap(keys, lookup)
}
# HH <- mat2hash(M)

# Lookup hashmap ----------------------------------------------------------
gen.key <- function(mu, nu) paste0(.f(mu), ',', .f(nu))
# (Also updating new key pairs).
hash.lookup <- function(mus, nus){
    keys <- gen.key(mus, nus)
    rtn <- HH[[keys]]
    w <- which(is.na(rtn))
    if(length(w)){ # Lazy --> Assumes nu is scalar value (i.e. one for all...)
        rtn[w] <- mapply(function(i) poly_solve(mus[i], nus), i = w, SIMPLIFY = T)
        HH$insert(keys[w], rtn[w])
    }
    rtn
}



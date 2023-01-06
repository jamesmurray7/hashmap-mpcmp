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
HH <- mat2hash(M)

# Lookup hashmap ----------------------------------------------------------
gen.key <- function(mu, nu) paste0(.f(mu), ',', .f(nu))

hash.lookup <- function(mu, nu){
    keys <- gen.key(mu, nu)
    HH[[keys]]
}

# Update hashmap ----------------------------------------------------------
# May be pointless?
check.current.hash <- function(H, mu, nu){

}

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

# Update hashmap ----------------------------------------------------------
# Maybe a little pointless? -->
# Check if current hashmap contains mus,nus.
check.current.hash <- function(keys) HH$has_keys(keys)
# Insert new values.
update.hash <- function(new.mus, new.nus){
    HH$insert(gen.key(new.mus, new.nus), c(poly_mat(new.mus, new.nus)))
    HH <<- HH
}

# Lookup hashmap ----------------------------------------------------------
gen.key <- function(mu, nu) paste0(.f(mu), ',', .f(nu))
hash.lookup <- function(mu, nu){
    keys <- gen.key(mu, nu)
    check <- check.current.hash(keys)
    if(!all(check)){
        keys
    }else{
        out <- HH[[keys]]
    }

    out
}



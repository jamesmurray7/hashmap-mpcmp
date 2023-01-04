# Trying out with hashmap -------------------------------------------------
mus <- nus <- seq(0.01, 2.00, 0.01)
M <- poly_mat(mus,nus)
# Dimnames
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

# Seems to work quite well!
HH[['1.52,1.20']] # poly_solve(1.52, 1.20)
HH[['0.01,0.01']] # poly_solve(0.01, 0.01)
HH[[c("1.00,1.01", "1.92,0.20", "0.20,1.24")]]
poly_solve(1,1.01); poly_solve(1.92,.2); poly_solve(.2,1.24)

# grid lookup vs. Benchmark -----------------------------------------------



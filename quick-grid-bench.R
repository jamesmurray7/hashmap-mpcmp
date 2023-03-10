# Proof of concept on 200; try full thing later.
mulong <- nulong <- seq(0.01, 10.00, 0.01)
microbenchmark::microbenchmark(
    `R` = {
        a <- matrix(0,200,200)
        for(i in 1:200){ # in actuality this uses foreach!
            for (j in 1:200){
                roots <- polyroot(((0:10) - mulong[i])/(factorial(0:10)^nulong[j]))
                lam <- Re(roots[Re(roots) >= 0 & zapsmall(Im(roots), 2) ==  0])
                # (Failure line here)
                # if (length(lam) == 0){lam <- grid_1K[i - 1, j]}
                a[i,j] <- lam
            }
        }
    },
    `C++` = {b <- poly_mat(mulong[1:200], nulong[1:200])},
    times = 5L
)

par(mfrow = c(2,1))
image(a, main = expression('Pete'~lambda), xlab = expression(mu), ylab = expression(nu))
image(b, main = expression('C++'~lambda), xlab = expression(mu), ylab = expression(nu))
max(abs(a-b))

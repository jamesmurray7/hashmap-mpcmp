library(microbenchmark)

mus <- runif(1000, 1.34, 2.50)
nus <- 1.5

mus.for.grid <- seq(min(mus) - 0.1, max(mus) + 0.1, .01) 
nus.for.grid <- seq(nus - .1, nus + .1, .01)

M <- generate.grid(mus.for.grid, nus.for.grid)
HH <- mat2hash(M)

# Function for lookups which are ALWAYS inside the grid...
benchN.inside <- function(Nmu, times = 1000, Nnu = 1L){
  mus <- runif(Nmu, min(mus.for.grid), max(mus.for.grid))
  nus <- runif(Nnu, min(nus.for.grid), max(nus.for.grid))
  
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
benchN(1e3, 1e3)
benchN(1e5, 1e2)

# Function for lookups which are CAN BE outside the grid...
benchN.outside <- function(Nmu, times = 1000, Nnu = 1L){
  mus <- runif(Nmu, min(mus.for.grid) - 0.2, max(mus.for.grid) + 0.2)
  nus <- runif(Nnu, min(nus.for.grid) - 0.2, max(nus.for.grid) + 0.2)
  
  .M <- M; .HH <- HH
  
  bench <- microbenchmark(
    `grid` = {
      G <- grid.lookup(mus, nus, update = T)
      M <<- .M
    },
    `hash` = {
      H <- hash.lookup(mus, nus)
      HH <<- .HH
    },
    times = times
  )
  bench
}

benchN.outside(1e3,100)

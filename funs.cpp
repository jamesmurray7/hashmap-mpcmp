// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;

// Generate polyroot grid -------------------------------------------------
// Function to solve
// [[Rcpp::export]]
vec poly_fn(const double mu, const double nu){
    vec a = reverse(regspace(0, 10));    // Opposite order to R's polyroot.
    return (a - mu)/(pow(tgamma(a + 1.), nu));
}

// [[Rcpp::export]]
double poly_solve(const double mu, const double nu){
    vec P = poly_fn(mu, nu);
    cx_vec C = roots(P);
    vec R = real(C), I = imag(C);
    uvec v = find(R >= 0. && I == 0.);
    return as_scalar(R.elem(v));
}

// Create grid (quite verbose!)
// [[Rcpp::export]]
mat poly_mat(const vec& mu, const vec& nu){
    int m = mu.size(), n = nu.size();
    Rcout << "Creating " << mu.size() << " x " << nu.size() << " matrix of lambda values..." << std::endl;
    // mu: row; nu: column;
    mat M = zeros<mat>(m, n);

    Rcout << std::endl;

    for(uword mm = 0; mm < m; mm++){
        for(uword nn = 0; nn < n; nn++){
            double lam = poly_solve(mu.at(mm), nu.at(nn));
            M.at(mm,nn) = lam;
        }
        if(mm % 100 == 0){
            Rcout << "mu: " << mm << "done." << std::endl;
        }
    }

    Rcout << std::endl << "Done!" << std::endl << std::endl;

    return M;
}


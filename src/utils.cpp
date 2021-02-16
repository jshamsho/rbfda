#include "utils.h"

arma::uvec arma_mod(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

arma::vec bayesreg(arma::vec const &b, arma::mat const &Q) {
  arma::vec z, w, v;
  z = arma::randn<arma::vec>(b.n_elem);
  arma::mat chol = arma::chol(Q, "lower");
  w = arma::solve(arma::trimatl(chol), b, arma::solve_opts::fast);
  v = arma::solve(arma::trimatu(chol.t()), w + z, arma::solve_opts::fast);
  return v;
}

arma::vec bayesreg_orth(arma::vec const &b, arma::mat const &Q, arma::mat &lin_constr) {
  arma::vec z, w, v;
  arma::mat chol, Cbar, Ctilde;
  z = arma::randn<arma::vec>(b.n_elem);
  chol = arma::chol(Q, "lower");
  w = arma::solve(arma::trimatl(chol), b, arma::solve_opts::fast);
  v = arma::solve(arma::trimatu(chol.t()), w + z, arma::solve_opts::fast);
  Cbar = arma::solve(arma::trimatl(chol), lin_constr.t(), arma::solve_opts::fast);
  Ctilde = arma::solve(arma::trimatu(chol.t()), Cbar, arma::solve_opts::fast);
  v = v - Ctilde * arma::solve(lin_constr * Ctilde, lin_constr * v);
  
  return v;
}

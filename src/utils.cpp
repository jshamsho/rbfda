#include "utils.h"

arma::uvec arma_mod(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

arma::vec bayesreg(arma::vec const &b, arma::mat const &Q) {
  arma::vec z, w, mu, v;
  arma::mat chol = arma::chol(Q, "lower");
  w = solve(arma::trimatl(chol), b, arma::solve_opts::fast);
  mu = solve(arma::trimatu(chol.t()), w, arma::solve_opts::fast);
  z = arma::randn<arma::vec>(b.n_elem);
  v = arma::solve(arma::trimatu(chol.t()), z, arma::solve_opts::fast);
  return mu + v;
}
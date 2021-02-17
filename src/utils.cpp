#include "utils.h"
static double const log2pi = std::log(2.0 * M_PI);

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

arma::vec bayesreg_orth(arma::vec const &b, arma::mat const &Q, arma::mat const &lin_constr) {
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


/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

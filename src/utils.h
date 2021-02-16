#ifndef MODSTRING_H1
#define MODSTRING_H1

#include <RcppArmadillo.h>
arma::uvec arma_mod(arma::uvec, arma::uword);
arma::vec bayesreg(arma::vec const &b, arma::mat const &Q);
arma::vec bayesreg_orth(arma::vec const &b, arma::mat const &Q, arma::mat &lin_constr);
#endif
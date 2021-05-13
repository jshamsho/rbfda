#ifndef MODSTRING_H1
#define MODSTRING_H1

#include <RcppArmadillo.h>
#include "data.h"
class Data;

arma::uvec arma_mod(arma::uvec, arma::uword);
arma::vec bayesreg(arma::vec const &b, arma::mat const &Q);
arma::vec bayesreg_orth(arma::vec const &b,
                        arma::mat const &Q,
                        arma::mat const &lin_constr);
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd);
double get_delta_eta_density(arma::mat& delta_eta1, arma::mat& delta_eta2,
                             arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                             arma::mat& design);
#endif
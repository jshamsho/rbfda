#ifndef MODSTRING_H1
#define MODSTRING_H1

#include <RcppArmadillo.h>
arma::uvec arma_mod(arma::uvec, arma::uword);
arma::vec bayesreg(arma::vec const &b, arma::mat const &Q);
arma::vec bayesreg_orth(arma::vec const &b,
                        arma::mat const &Q,
                        arma::mat const &lin_constr);
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd);
#endif
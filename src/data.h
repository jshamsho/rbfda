#ifndef DATA_H
#define DATA_H
#include <RcppArmadillo.h>
#include "utils.h"

class Data {
public:
  arma::uword basisdim, ldim, nsub, nreg, nt, designdim;
  arma::uvec missing, missing_sub, missing_time, missing_reg;
  arma::mat response;
  arma::mat design;
  arma::mat basis;
  arma::mat penalty;
  arma::uword iter, thin, burnin;
  arma::vec time;
  Data() {};
  Data(arma::mat&, arma::mat&, 
       arma::mat&,
       arma::vec&, arma::mat&,
       arma::uword, arma::uword, 
       arma::uword, arma::uword);
  Rcpp::List write_data();
};
//

#endif
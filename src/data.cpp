#include "data.h"

Data::Data(arma::mat& response, arma::mat& design,
           arma::mat& basis, arma::vec& time, 
           arma::field<arma::mat>& penalties,
           arma::uvec& indices, 
           arma::uword kdim, arma::uword iter,
           arma::uword burnin, arma::uword thin)
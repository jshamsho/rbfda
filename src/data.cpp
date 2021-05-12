#include "data.h"

Data::Data(arma::mat& response, arma::mat& design, 
           arma::mat& basis,
           arma::vec& time, arma::mat& penalty,
           arma::uword ldim, arma::uword iter,
           arma::uword burnin, arma::uword thin) {
  this->response = response;
  this->design = design;
  this->designdim = design.n_cols;
  this->basis= basis;
  this->penalty = penalty;
  this->penalty_rank = arma::rank(penalty);
  this->ldim = ldim;
  this->basisdim = basis.n_cols;
  this->time = time;
  this->burnin = burnin;
  this->iter = iter;
  this->thin = thin;
  this->nreg = response.n_cols;
  this->nt = time.n_elem;
  this->nsub = design.n_rows;
  this->missing = arma::find_nonfinite(this->response);
  this->missing_sub = arma_mod(this->missing, this->nt * this->nsub);
  this->missing_reg = arma::floor(this->missing / (this->nt * this->nsub));
  this->c = std::min(nreg, ldim);
}

Rcpp::List Data::write_data() {
  return Rcpp::List::create(Rcpp::Named("response", response),
                            Rcpp::Named("basis", basis),
                            Rcpp::Named("time", time),
                            Rcpp::Named("ldim", ldim),
                            Rcpp::Named("missing_subjects", missing_sub),
                            Rcpp::Named("missing_time", missing_time),
                            Rcpp::Named("design_mean", design),
                            Rcpp::Named("nsub", nsub),
                            Rcpp::Named("nreg", nreg));
}
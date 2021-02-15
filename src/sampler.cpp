#include "sampler.h"

// [[Rcpp::depends(RcppArmadillo)]]
/*
Sampler::Sampler(Data& dat, Parameters& pars, Transformations& transf) {
  this->dat = dat;
  this->pars = pars;
  this->transf = transf;
}
*/
Rcpp::List Sampler::write_data() {
  return(dat.write_data());
}
/*
Rcpp::List Sampler::write_control() {
  return(Rcpp::List::create(
      Rcpp::Named("iterations", dat.iter),
      Rcpp::Named("thin", dat.thin),
      Rcpp::Named("burnin", dat.burnin)
  ));
}

void Sampler::sample() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort()) {
        Rcpp::Rcout << "MCMC terminated by user\n";
        goto stop;
      }
      
    }
    progress_bar.increment();
    write_samples();
  }
  stop:
    NULL;
}
*/

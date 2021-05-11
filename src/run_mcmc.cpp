#include <RcppArmadillo.h>
#include <iostream>
#include "data.h"
#include "parameters.h"
#include "transformations.h"
#include "sampler.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Run Markov-Chain Monte-Carlo
//'
//' Generate samples from the posterior distribution using Gibbs sampling and 
//' Metropolis-Hastings steps
//' @param response N x T response matrix, where N is number of subjects and T is number
//' of time points. Values can be NA if there's missing data
//' @param design N x d_{1} design matrix for mean structure
//' @param basis T x p basis matrix for functional domain
//' @param time vector of time points
//' @param penalties List of smoothing penalties for mean structure
//' @param indices_mean Maps penalties in the mean structure to beta coefficients
//' @param kdim Dimension of latent subspace
//' @param iter Number of iterations to run
//' @param burnin Number of iterations to use as burn-in. This is only relevant when
//' passing the returned object into post-processing functions
//' @param thin Thinning defaulting to 1
//' @export run_mcmc
//' @return A List containing 3 lists including data, control, and samples.
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat response, arma::mat design,
                    arma::mat basis, arma::vec time,
                    arma::mat penalty, 
                    arma::uword ldim, arma::uword iter, arma::uword burnin,
                    arma::uword thin = 1, Rcpp::Nullable<Rcpp::List> init_ = R_NilValue,
                    std::string covstruct = "partial") {
  
  Data dat(response, design,
           basis, time,
           penalty, ldim,
           iter, burnin, thin);
  Sampler* mysampler = SamplerFactory::new_mcmc(covstruct, dat, init_);
  mysampler->sample();
  Rcpp::List return_me;
  return_me["data"] = mysampler->write_data();
  return_me["samples"] = mysampler->get_samples();
  return_me["control"] = mysampler->write_control();
  return return_me;
}
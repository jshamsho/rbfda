#ifndef SAMPLER_H
#define SAMPLER_H

#include <RcppArmadillo.h>
#include "data.h"
#include "transformations.h"
#include "parameters.h"
#include <progress.hpp>
class Data;
class Transformations;
class Parameters;
class Sampler
{
  public:
    Data dat;
    Parameters pars;
    Transformations transf;
    Sampler(Data& dat, Parameters& pars, Transformations& transf) :
      dat(dat), pars(pars), transf(transf) {}
    // Sampler(Data&, Parameters& pars, Transformations& transf);
    void sample();
    void write_samples();
    Rcpp::List get_samples();
    Rcpp::List write_data();
    Rcpp::List write_control();
};

#endif

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
  Parameters *pars;
  Transformations transf;
  // parameterized constructor
  Sampler(Data& dat, Parameters& pars, Transformations& transf) :
    dat(dat), pars(&pars), transf(transf) {}
  arma::uword current_iter = 0;
  virtual void sample() = 0;
  virtual void write_samples() = 0;
  virtual Rcpp::List get_samples() = 0;
  Rcpp::List write_data();
  Rcpp::List write_control();
  ~Sampler(){}
};

class SamplerPartial : public Sampler
{
  public:
    SamplerPartial(Data& dat, Parameters& pars, Transformations& transf) :
    Sampler(dat, pars, transf) {}
    void sample();
    void write_samples();
    Rcpp::List get_samples();
    ~SamplerPartial() {}
};

class SamplerFactory {
public:
  static Sampler *new_mcmc(std::string covstruct, Data& dat, Parameters& pars, Transformations& transf) {
    if(covstruct == "partial") return new SamplerPartial(dat, pars, transf);
    // if(type == "unequal") return new SamplerUnequal(dat, pars, transf);
    return nullptr;
  }
};

#endif

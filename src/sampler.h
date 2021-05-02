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
  // Parameters *pars;
  // Transformations transf;
  // parameterized constructor
  // Sampler(Data& dat, Parameters& pars, Transformations& transf) :
    // dat(dat), pars(&pars), transf(transf) {}
  Sampler() {};
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
    ParametersPartial pars;
    TransformationsPartial transf;
    SamplerPartial(Data&, Rcpp::Nullable<Rcpp::List>);
    void sample();
    void write_samples();
    Rcpp::List get_samples();
    ~SamplerPartial() {}
};

class SamplerWeak : public Sampler
{
public:
  ParametersWeak pars;
  Transformations transf;
  SamplerWeak(Data& dat, Rcpp::Nullable<Rcpp::List>);
  void sample();
  void write_samples();
  Rcpp::List get_samples();
  ~SamplerWeak() {}
};
class SamplerFactory {
public:
  static Sampler *new_mcmc(std::string covstruct, Data& dat, Rcpp::Nullable<Rcpp::List> init_) {
    if(covstruct == "partial") {
      return new SamplerPartial(dat, init_); 
    } else if (covstruct == "weak") {
      return new SamplerWeak(dat, init_);
    }
    return nullptr;
  }
};

#endif

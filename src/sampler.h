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
  TransformationsWeak transf;
  SamplerWeak(Data& dat, Rcpp::Nullable<Rcpp::List>);
  void sample();
  void write_samples();
  arma::vec get_grad(arma::vec&);
  arma::vec p, q;
  arma::mat leapfrog(arma::uword num_steps, double step_size);
  void update_delta_eta(arma::uword num_steps, double step_size);
  double get_density(arma::vec& p);
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

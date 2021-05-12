#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "data.h"
#include "transformations.h"

class Transformations;
class TransformationsPartial;
class TransformationsWeak;

class Parameters
{
  public:
    // psi: Loading matrix for functional domain
    // phi: Loading matrix for region domain
    // sigmasqetai: region by # of latent factor random component precision
    // eta: nreg * nsub by ldim subject scores
    // delta: Subject-region-time specific precisions, size nt * nsub by nreg
    //   see https://jrnold.github.io/bayesian_notes/robust-regression.html
    // omega: Region-specific precisions
    // zeta: Smoothing parameters for psi with uniform prior on SD scale
    //   Corresponds to Gamma(-.5, 0) on precision scale
    // nup: Region specific degrees of freedom=
    // 
    double prior_zeta_shape = -.5;
    double prior_zeta_rate = 0;
    double prior_omega_shape = 0;
    double prior_omega_rate = 0;
    double posterior_omega_shape;
    double delta_beta_nu = 4;
    double rho_shape1 = 1, rho_shape2 = 1;
    double old_logpost = -arma::datum::inf;
    double rho, alpha;
    double nu;
    arma::vec omega, zeta, rho_container, alpha_container,
      nu_container, delta_lambda;
    arma::mat lambda, sigmasqeta, sigmasqetai, delta, eta, xi_eta,
      beta, delta_beta, phi0, tau_phi0, xi_lambda, delta_phi;
    arma::cube phi, xi_lambda_container;
    arma::mat omega_container, zeta_container;
    arma::cube lambda_container, eta_container, sigmasqetai_container, 
      sigmasqeta_container, xi_eta_container, beta_container, delta_beta_container,
      phi0_container, tau_phi0_container, delta_lambda_container, 
      xi_phi, delta_phi_container;
    arma::field<arma::cube> xi_phi_container;
    
    // omega_container = arma::mat(dat.nreg, dat.iter);
    // zeta_container = arma::mat(dat.ldim, dat.iter);
    // delta_eta_container = arma::mat(dat.nreg, dat.ldim, dat.iter);
    double prior_shape = 0;
    Parameters(const Data&, Rcpp::Nullable<Rcpp::List>);
    Parameters() {};
    // virtual void update_lambda(const Data&, Transformations&) = 0;
    // virtual void update_phi(const Data&, Transformations&) = 0;
    // virtual void update_eta(const Data&, Transformations&) = 0;
    // virtual void update_omega(const Data&, Transformations&) = 0;
    
    void update_delta(const Data&, Transformations&);
    void update_xi_eta(const Data&, Transformations&);
    // void update_delta_eta(const Data&, Transformations&);
    void update_beta(const Data&, Transformations&);
    void update_delta_beta(const Data&, Transformations&);
    void update_zeta(const Data&, Transformations&);
    void update_nu(const Data&, Transformations);
    // void update_a123(const Data&, Transformations&);
    ~Parameters(){}
};

class ParametersPartial : public Parameters
{
public:
  double a1, a2, a3;
  arma::vec a1_container, a2_container, a3_container;
  arma::cube phi;
  arma::field<arma::cube> phi_container;
  arma::mat delta_eta;
  arma::cube delta_eta_container;
  ParametersPartial() {}
  ParametersPartial(const Data& dat, Rcpp::Nullable<Rcpp::List> init_);
  void update_eta(const Data&, Transformations&);
  void update_delta_eta(const Data&, Transformations&);
  void update_a123(const Data&);
  void update_lambda(const Data&, Transformations&);
  void update_phi(const Data&, TransformationsPartial&);
  void update_omega(const Data&, TransformationsPartial&);
  ~ParametersPartial() {}
};

class ParametersWeak : public Parameters
{
public:
  double a1, a2, a3, a4;
  arma::mat phi;
  arma::cube phi_container;
  arma::vec delta_eta1, delta_eta2, a1_container, a2_container, a3_container,
    a4_container;
  arma::mat delta_eta1_container, delta_eta2_container;
  ParametersWeak() {}
  ParametersWeak(const Data& dat, Rcpp::Nullable<Rcpp::List> init_);
  void update_eta(const Data&, TransformationsWeak&);
  void update_delta_eta1(const Data&, TransformationsWeak&);
  void update_delta_eta2(const Data&, TransformationsWeak&);
  // void update_delta_eta(const Data&, Transformations&);
  void update_a123(const Data&);
  void update_lambda(const Data&, TransformationsWeak&);
  void update_phi(const Data&, TransformationsWeak&);
  void update_omega(const Data&, TransformationsWeak&);
  ~ParametersWeak() {}
};

class ParametersFactory {
public:
  static Parameters *new_pars(std::string covstruct, Data& dat, Rcpp::Nullable<Rcpp::List> init_) {
    if(covstruct == "partial") return new ParametersPartial(dat, init_);
    if(covstruct == "weak") return new ParametersWeak(dat, init_);
    // if(type == "unequal") return new SamplerUnequal(dat, pars, transf);
    return nullptr;
  }
};
#endif
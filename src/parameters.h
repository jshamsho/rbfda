#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "data.h"
#include "transformations.h"

class Transformations;

class Parameters {
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
    double delta_nu = 5;
    double delta_eta_nu = 5;
    double delta_beta_nu = 4;
    arma::vec omega, nup, zeta;
    arma::mat lambda, sigmasqeta, sigmasqetai, delta, eta, xi_eta, delta_eta,
      beta, delta_beta;
    arma::cube phi;
    arma::mat omega_container, zeta_container;
    arma::cube lambda_container, eta_container, sigmasqetai_container, 
      sigmasqeta_container, xi_eta_container, beta_container, delta_beta_container,
      delta_eta_container;
    arma::field<arma::cube> phi_container;
    
    // omega_container = arma::mat(dat.nreg, dat.iter);
    // zeta_container = arma::mat(dat.ldim, dat.iter);
    // delta_eta_container = arma::mat(dat.nreg, dat.ldim, dat.iter);
    double prior_shape = 0;
    Parameters(Data&);
    Parameters() {};
    void update_omega(Data&, Transformations&);
    void update_delta(Data&, Transformations&);
    void update_eta(Data&, Transformations&);
    void update_xi_eta(Data&, Transformations&);
    void update_delta_eta(Data&, Transformations&);
    void update_beta(const Data&, Transformations&);
    void update_delta_beta(const Data&, Transformations&);
    void update_lambda(const Data&, Transformations&);
    void update_zeta(const Data&, Transformations&);
    void update_phi(const Data&, Transformations&);
};
#endif
#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <RcppArmadillo.h>
#include "parameters.h"
#include "data.h"
//#include "GaussTransformations.h"

class Parameters;

class Transformations {
public:
  Transformations() {};
  Transformations(Data&, Parameters&);
  void initialize_fit(Data&, Parameters&);
  void complete_response(Data&, Parameters&);
  arma::mat fit, fit_latent, fit_eta, psi, btb, bty, psi_lin_constr,
  phi_lin_constr, C_rho, ones_mat, delta_eta_cumprod, delta_eta_cumprod_init;
  
};

#endif

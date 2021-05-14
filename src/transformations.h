#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <RcppArmadillo.h>
#include "parameters.h"
#include "data.h"
//#include "GaussTransformations.h"

class Parameters;
class ParametersPartial;
class ParametersWeak;

class Transformations {
public:
  Transformations() {};
  Transformations(Data&, Parameters&);
  // virtual void initialize_fit(Data&, Parameters&) = 0;
  void complete_response(Data&, Parameters&);
  arma::mat fit, fit_latent, fit_eta, psi, btb, bty, psi_lin_constr,
    phi_lin_constr, C_rho, ones_mat, delta_eta_cumprod_init;
};

class TransformationsPartial : public Transformations {
public:
  void initialize_fit(Data&, ParametersPartial&);
  TransformationsPartial() {}
  TransformationsPartial(Data&, ParametersPartial&);
};

class TransformationsWeak : public Transformations {
public:
  void initialize_fit(Data&, ParametersWeak&);
  TransformationsWeak() {}
  TransformationsWeak(Data&, ParametersWeak&);
};
#endif

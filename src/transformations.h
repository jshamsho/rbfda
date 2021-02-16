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
  arma::mat fit, fit_latent, fit_eta, psi, btb, bty;
  
};

#endif

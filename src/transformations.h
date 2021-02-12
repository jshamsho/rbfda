#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <RcppArmadillo.h>
#include "parameters.h"
#include "data.h"
//#include "GaussTransformations.h"

class Parameters;

class Transformations {
public:
  Transformations(Data&, Parameters&);
  arma::mat fit;
  
};

#endif

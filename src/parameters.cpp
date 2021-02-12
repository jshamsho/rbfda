#include "parameters.h"

Parameters::Parameters(Data& dat) {
  omega = arma::vec(dat.nreg);
  posterior_omega_shape = dat.nt * dat.nsub / 2. + prior_shape;
  delta = arma::mat(dat.nt * dat.nsub, dat.nreg, arma::fill::ones);
  psi = arma::mat(dat.nt, dat.ldim, arma::fill::randn);
  phi = arma::cube(dat.nreg, dat.nreg, dat.ldim, arma::fill::zeros);
  eta = arma::mat(dat.nreg * dat.nsub, dat.ldim, arma::fill::randn);
  sigmasqetai = arma::mat(dat.nreg, dat.ldim);
}

void Parameters::update_omega(Data& dat, Transformations& transf) {
  double rate;
  for (arma::uword r = 0; r < dat.nreg; r++) {
      rate = .5 * arma::dot(delta.col(r), arma::square(dat.response.col(r) -
        transf.fit.col(r)));
    omega(r) = R::rgamma(posterior_omega_shape, 1. / (prior_omega_rate + rate));
  }
}

void Parameters::update_delta(Data& dat, Transformations& transf) {
  double shape = .5 + delta_nu;
  double rate;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword i = 0; i < dat.nt * dat.nsub; i++) {
      rate = delta_nu + .5 * std::pow(dat.response(i, r) - transf.fit(i, r), 2) *
        omega(r);
      delta(i, r) = R::rgamma(shape, 1. / rate);
    }
  }
}


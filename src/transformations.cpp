#include "transformations.h"

Transformations::Transformations(Data& dat, Parameters& pars) {
  btb = dat.basis.t() * dat.basis;
  fit = arma::mat(dat.nsub * dat.nt, dat.nreg, arma::fill::zeros);
  fit_eta = arma::mat(dat.nreg, dat.ldim, arma::fill::randn);
  psi = arma::mat(dat.nt, dat.ldim, arma::fill::randn);
  initialize_fit(dat, pars);
}

void Transformations::initialize_fit(Data& dat, Parameters& pars) {
  arma::uword first = 0;
  arma::uword last = dat.nt - 1;
  arma::uword first_eta = 0;
  arma::uword last_eta = dat.nreg - 1;
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      fit.rows(first, last) =
        fit.rows(first, last) + psi.col(l) * 
        (arma::mat(pars.eta.rows(first_eta, last_eta)).col(l).t() * pars.phi.slice(l).t());
      first = i * dat.nt;
      last = (i + 1) * dat.nt - 1;
      first_eta = i * dat.nreg;
      last_eta = (i + 1) * dat.nreg - 1;
    }
  }
}


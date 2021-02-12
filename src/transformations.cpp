#include "transformations.h"

Transformations::Transformations(Data& dat, Parameters& pars) {
  fit = arma::mat(dat.nsub * dat.nt, dat.nreg, arma::fill::zeros);
  arma::uword first = 0;
  arma::uword last = dat.nt - 1;
  arma::uword first_eta = 0;
  arma::uword last_eta = dat.nreg - 1;
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      fit.rows(first, last) =
        fit.rows(first, last) + pars.psi.col(l) * 
        (arma::mat(pars.eta.rows(first_eta, last_eta)).col(l).t() * pars.phi.slice(l).t());
            // (pars.phi.slice(l) * pars.eta.slice(i).col(l)) * pars.psi.col(l).t());
      first = i * dat.nt;
      last = (i + 1) * dat.nt - 1;
      first_eta = i * dat.nreg;
      last_eta = (i + 1) * dat.nreg - 1;
    }
  }
}
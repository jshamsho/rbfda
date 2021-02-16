#include "transformations.h"

Transformations::Transformations(Data& dat, Parameters& pars) {
  btb = dat.basis.t() * dat.basis;
  bty = arma::mat(dat.basisdim * dat.nsub, dat.nreg);
  fit = arma::mat(dat.nsub * dat.nt, dat.nreg, arma::fill::zeros);
  fit_latent = arma::mat(dat.nsub * dat.basisdim, dat.nreg, arma::fill::zeros);
  fit_eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::randn);
  psi = dat.basis * pars.lambda;
  lin_constr = arma::mat(dat.ldim, dat.basisdim);
  initialize_fit(dat, pars);
}

void Transformations::initialize_fit(Data& dat, Parameters& pars) {
  arma::uword first = 0;
  arma::uword last = dat.nt - 1;
  arma::uword first_eta = 0;
  arma::uword last_eta = dat.nreg - 1;
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      first = i * dat.nt;
      last = (i + 1) * dat.nt - 1;
      first_eta = i * dat.nreg;
      last_eta = (i + 1) * dat.nreg - 1;
      fit.rows(first, last) =
        fit.rows(first, last) + psi.col(l) * 
        (arma::mat(pars.eta.rows(first_eta, last_eta)).col(l).t() * pars.phi.slice(l).t());
    }
  }
  complete_response(dat, pars);
  for (arma::uword i = 0; i < dat.nsub; i++) {
    bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) = dat.basis.t() * dat.response.rows(i * dat.nt, (i + 1) * dat.nt - 1);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      fit_latent.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) = 
        fit_latent.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) + 
        pars.lambda.col(l) * pars.eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() * pars.phi.slice(l).t();
    }
  }
  
  for (arma::uword l = 0; l < dat.ldim; l++) {
    lin_constr.row(l) = pars.lambda.col(l).t() * btb;
  }
}

void Transformations::complete_response(Data& dat, Parameters& pars) {
  arma::uword r;
  for (arma::uword z = 0; z < dat.missing.n_elem; z++) {
    r = dat.missing_reg(z);
    dat.response(dat.missing(z)) = fit(dat.missing(z)) + R::rnorm(0, std::pow(pars.omega(r), -.5));
  }
  /*
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      for (arma::uword t = 0; t < dat.nt; t++) {
        if (dat.missing(i * dat.nt + ))
        dat.response.row(i * dat.nt + t).col(r) = fit.row(i * dat.nt + t).col(r) + 
          R::rnorm(0, std::pow(pars.omega(r), -.5));
      } 
    }
  }*/
}
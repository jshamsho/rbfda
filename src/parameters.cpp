#include "parameters.h"

Parameters::Parameters(Data& dat) {
  omega = arma::vec(dat.nreg, arma::fill::ones);
  posterior_omega_shape = dat.nt * dat.nsub / 2. + prior_shape;
  delta = arma::mat(dat.nt * dat.nsub, dat.nreg, arma::fill::ones);
  zeta = arma::vec(dat.ldim, arma::fill::ones);
  lambda = arma::mat(dat.basisdim, dat.ldim, arma::fill::randn);
  phi = arma::cube(dat.nreg, dat.nreg, dat.ldim, arma::fill::randn);
  eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::randn);
  sigmasqetai = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
  sigmasqeta = arma::mat(dat.nreg, dat.ldim, arma::fill::ones);
  xi_eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
  beta = arma::mat(dat.designdim * dat.nreg, dat.ldim, arma::fill::randn);
  delta_beta = arma::mat(dat.nreg * dat.designdim, dat.ldim, arma::fill::ones);
  delta_eta = arma::mat(dat.nreg, dat.ldim, arma::fill::ones);
}

void Parameters::update_omega(Data& dat, Transformations& transf) {
  double rate;
  for (arma::uword r = 0; r < dat.nreg; r++) {
      rate = .5 * arma::dot(delta.col(r), arma::square(dat.response.col(r) -
        transf.fit.col(r)));
    omega(r) = R::rgamma(posterior_omega_shape, 1. / (prior_omega_rate + rate));
  }
}

/*
 * Working with normal likelihood for now. Leads to simpler working likelihood
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
*/

// eta_il = beta * xi + eps
// (nreg x 1), (nreg x d), (d x 1), (nreg x 1)
void Parameters::update_eta(Data& dat, Transformations& transf) {
  arma::uword first, last, first_eta, last_eta;
  arma::vec b, yt;
  arma::mat Q, beta_mat, diagsigma, diagomega = arma::diagmat(omega);
  for (arma::uword i = 0; i < dat.nsub; i++) {
    first = i * dat.nt, last = (i + 1) * dat.nt - 1;
    first_eta = i * dat.nreg, last_eta = i * dat.nreg + dat.nreg - 1;
    for (arma::uword l = 0; l < dat.ldim; l++) {
      diagsigma = arma::diagmat(sigmasqetai.rows(first_eta, last_eta).col(l));
      yt = arma::trans((transf.psi.col(l).t() * dat.response.rows(first, last)) * phi.slice(l));
      beta_mat = arma::reshape(beta.col(l), dat.designdim, dat.nreg).t();
      b = phi.slice(l).t() * diagomega * phi.slice(l) * yt + 
        diagsigma * (beta_mat * dat.design.row(i).t());
      Q = phi.slice(l).t() * diagomega * phi.slice(l) + diagsigma;
      eta.rows(first_eta, last_eta).col(l) = bayesreg(b, Q);
    }
  }
}

void Parameters::update_xi_eta(Data& dat, Transformations& transf) {
  double shape = .5 + delta_eta_nu;
  double rate;
  arma::uword idx;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      for (arma::uword i = 0; i < dat.nsub; i++) {
        idx = i * dat.nreg + r;
        rate = delta_eta_nu + std::pow(eta(idx, l) - transf.fit_eta(idx, l), 2) * sigmasqeta(r, l);
        xi_eta(i * dat.nreg + r, l) = R::rgamma(shape, 1. / rate);
      }
    }
  }
}

void Parameters::update_delta_eta(Data& dat, Transformations& transf) {/*
  double rate;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      rate = 
    }
    rate = .5 * arma::dot(delta.col(r), arma::square(dat.response.col(r) -
      transf.fit.col(r)));
    omega(r) = R::rgamma(posterior_omega_shape, 1. / (prior_omega_rate + rate));
  }*/
}

void Parameters::update_beta(const Data& dat, Transformations& transf) {
  arma::uword first, last;
  arma::uvec r_ind;
  arma::vec eta_vec, b;
  arma::mat Q, db, sig_eta;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      first = r * dat.designdim;
      last = (r + 1) * dat.designdim - 1;
      r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
      eta_vec = arma::vec(eta.col(l)).elem(r_ind);
      sig_eta = arma::diagmat(arma::vec(sigmasqetai.col(l)).elem(r_ind));
      db = arma::diagmat(delta_beta.col(l).rows(first, last));
      Q = dat.design.t() * sig_eta * dat.design + db;
      b = dat.design.t() * (sig_eta * eta_vec);
      beta.col(l).rows(first, last) = bayesreg(b, Q);
    }
  }
}

void Parameters::update_delta_beta(const Data& dat, Transformations& transf) {
  double shape = delta_beta_nu + .5;
  double rate = delta_beta_nu;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword j = 0; j < dat.designdim * dat.nreg; j++) {
      rate = rate + .5 * std::pow(beta(j, l), 2.);
      delta_beta(j, l) = R::rgamma(shape, 1. / rate);
    }
  }
}

void Parameters::update_lambda(const Data& dat, Transformations& transf) {
  arma::mat eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  arma::vec eta_temp;
  arma::vec b = arma::zeros(dat.basisdim);
  arma::mat Q;
  arma::mat diagomega = arma::diagmat(omega);
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_sum = eta_sum + eta_temp * eta_temp.t();
      b = b + dat.basis.t() * dat.response.rows(i * dat.nt, (i + 1) * dat.nt - 1) * (diagomega * (phi.slice(l) * eta_temp));
      for (arma::uword ll = 0; ll < dat.ldim; ll++) {
        if (ll != l) {
          b = b - (dat.basis.t() * transf.psi.col(l)) * (eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() *
            (phi.slice(ll).t() * (diagomega * (phi.slice(l) * eta_temp))));
        }
      }
    }
    Q = arma::trace(phi.slice(l).t() * diagomega * phi.slice(l) * eta_sum) * transf.btb + zeta(l) * dat.penalty;
    lambda.col(l) = bayesreg(b, Q);
    b.zeros();
  }

}
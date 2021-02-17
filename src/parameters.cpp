#include "parameters.h"

Parameters::Parameters(Data& dat) {
  omega_container = arma::mat(dat.nreg, dat.iter);
  zeta_container = arma::mat(dat.ldim, dat.iter);
  lambda_container = arma::cube(dat.basisdim, dat.ldim, dat.iter);
  phi_container = arma::field<arma::cube>(dat.iter);
  eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqetai_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqeta_container = arma::cube(dat.nreg, dat.ldim, dat.iter);
  xi_eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  beta_container = arma::cube(dat.designdim * dat.nreg, dat.ldim, dat.iter);
  delta_beta_container = arma::cube(dat.nreg * dat.designdim, dat.ldim, dat.iter);
  delta_eta_container = arma::cube(dat.nreg, dat.ldim, dat.iter);
  rho_container = arma::vec(dat.iter);
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
      rate = prior_omega_rate +
        .5 * arma::dot(delta.col(r), arma::square(dat.response.col(r) -
        transf.fit.col(r)));
    omega(r) = R::rgamma(posterior_omega_shape, 1. / rate);
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

void Parameters::update_eta(Data& dat, Transformations& transf) {
  arma::uword first, last, first_eta, last_eta;
  arma::vec b, yt;
  arma::mat Q, beta_mat, diagsigma, diagomega = arma::diagmat(omega);
  transf.fit.zeros();
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
      transf.fit.rows(first, last) = transf.fit.rows(first, last) +
        transf.psi.col(l) * (eta.rows(first_eta, last_eta).col(l).t() * phi.slice(l).t());
    }
  }
}

void Parameters::update_xi_eta(Data& dat, Transformations& transf) {
  double shape = .5 + delta_eta_nu / 2.;
  double rate = 0;
  arma::uword idx;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      for (arma::uword i = 0; i < dat.nsub; i++) {
        idx = i * dat.nreg + r;
        rate = delta_eta_nu +
          std::pow(eta(idx, l) - transf.fit_eta(idx, l), 2) * sigmasqeta(r, l);
        xi_eta(idx, l) = R::rgamma(shape, 1. / rate);
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
      arma::vec(transf.fit_eta.col(l)).rows(r_ind) = dat.design * beta.col(l).rows(first, last);
    }
  }
}

void Parameters::update_delta_beta(const Data& dat, Transformations& transf) {
  double shape = delta_beta_nu / 2. + .5;
  double rate = delta_beta_nu / 2.;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword j = 0; j < dat.designdim * dat.nreg; j++) {
      rate = rate + .5 * std::pow(beta(j, l), 2.);
      delta_beta(j, l) = R::rgamma(shape, 1. / rate);
    }
  }
}

void Parameters::update_lambda(const Data& dat, Transformations& transf) {
  double psi_norm;
  arma::mat Q, eta_sum, eta_phi, diagomega;
  arma::vec b, eta_temp;
  eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  eta_phi = arma::mat(dat.nsub, dat.nreg);
  diagomega = arma::diagmat(omega);
  b = arma::vec(dat.basisdim);

  for (arma::uword l = 0; l < dat.ldim; l++) {
    b.zeros();
    eta_sum.zeros();

    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_sum = eta_sum + eta_temp * eta_temp.t();
      eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) -
        transf.psi.col(l) * eta_phi.row(i);
      b = b + (transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) -
        dat.basis.t() * transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1)) *
        diagomega * (eta_phi.row(i)).t();
    }
    Q = arma::trace(phi.slice(l).t() * diagomega * phi.slice(l) * eta_sum) * transf.btb + zeta(l) * dat.penalty;
    transf.psi_lin_constr.shed_row(l);
    lambda.col(l) = bayesreg_orth(b, Q, transf.psi_lin_constr);
    psi_norm = arma::norm(dat.basis * lambda.col(l));
    lambda.col(l) = lambda.col(l) / psi_norm;
    transf.psi.col(l) = dat.basis * lambda.col(l);
    eta.col(l) = eta.col(l) * psi_norm;
    transf.psi_lin_constr.insert_rows(l, transf.psi.col(l).t() * dat.basis);

    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
        transf.psi.col(l) * eta_phi.row(i);
    }
  }
}

void Parameters::update_zeta(const Data& dat, Transformations& transf) {
  double shape = prior_zeta_shape + .5 * dat.basisdim;
  double rate;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    rate = prior_zeta_rate + .5 *
      arma::as_scalar(lambda.col(l).t() * dat.penalty * lambda.col(l));
    zeta(l) = R::rgamma(shape, 1. / rate);
  }
}

void Parameters::update_phi(const Data& dat, Transformations& transf) {
  double norm;
  arma::uword first, last, idx;
  arma::vec b, phi_temp;
  arma::uvec r_ind;
  arma::mat Q, diagomega, diageta, eta_sum;
  b = arma::zeros(dat.nreg * dat.ldim);
  diagomega = arma::diagmat(omega);
  eta_sum = arma::zeros(dat.ldim, dat.ldim);
  transf.fit.zeros();
  arma::uvec constr_indices;// = arma::uvec(0);
  // Rcpp::Rcout << "head of func: \n" << transf.phi_lin_constr << "\n";
  for (arma::uword r = 0; r < dat.nreg; r++) {
    b.zeros();
    eta_sum.zeros();
    r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nreg - 1) * dat.ldim + r);
    for (arma::uword i = 0; i < dat.nsub; i++) {
      first = i * dat.nt, last = (i + 1) * dat.nt - 1;
      idx = i * dat.nreg + r;
      diageta = arma::diagmat(eta.row(idx));
      b = b + arma::vectorise(diagomega * dat.response.rows(first, last).t() *
        transf.psi * diageta);
      eta_sum = eta_sum + diageta * diageta;
    }
    Q = arma::kron(eta_sum, diagomega) +
      arma::kron(transf.C_rho, arma::eye(dat.nreg, dat.nreg));
    // transf.phi_lin_constr.shed_rows(r_ind);
    if (r > 0) {
      constr_indices =
        arma::sort(arma::join_cols(
          constr_indices,
          arma::regspace<arma::uvec>(
            r - 1, dat.nreg, (dat.ldim - 1) * dat.nreg + r)));
    }
    phi_temp = bayesreg_orth(b, Q, transf.phi_lin_constr.rows(constr_indices));
    

    for (arma::uword l = 0; l < dat.ldim; l++) {
      norm = arma::norm(phi_temp.rows(l * dat.nreg, (l + 1) * dat.nreg - 1));
      phi.slice(l).col(r) = phi_temp.rows(l * dat.nreg, (l + 1) * dat.nreg - 1) / norm;
      transf.phi_lin_constr.row(l * dat.nreg + r).cols(l * dat.nreg, (l + 1) * dat.nreg - 1) =
        phi.slice(l).col(r).t();
      for (arma::uword i = 0; i < dat.nsub; i++) {
        eta.row(i * dat.nreg + r).col(l) = eta.row(i * dat.nreg + r).col(l) * norm;
      }
    }
  }

  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
        transf.psi.col(l) * (eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() *
        phi.slice(l).t());
    }
  }
}

void Parameters::update_rho(const Data &dat, Transformations &transf) {
  double loglik = 0;
  arma::mat C_rho = rho * transf.ones_mat + (1 - rho) * arma::eye(dat.ldim, dat.ldim);
  arma::mat x = arma::mat(dat.nreg * dat.nreg, dat.ldim);
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword rr = 0; rr < dat.nreg; rr++) {
      x.row(rr + r * dat.nreg) = arma::vec(phi.tube(r, rr)).t();
    }
  }
  loglik = arma::accu(dmvnrm_arma_fast(x, arma::zeros<arma::rowvec>(dat.ldim), C_rho, true));
  // Rcpp::Rcout << "loglik: " << loglik << "\n";
}

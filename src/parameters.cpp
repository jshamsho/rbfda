#include "parameters.h"

ParametersPartial::ParametersPartial(const Data& dat, Rcpp::Nullable<Rcpp::List> init_ = R_NilValue) {
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
  a1_container = arma::vec(dat.iter);
  a2_container = arma::vec(dat.iter);
  a3_container = arma::vec(dat.iter);
  alpha_container = arma::vec(dat.iter);
  nu_container = arma::vec(dat.iter);
  phi = arma::cube(dat.nreg, dat.nreg, dat.ldim, arma::fill::zeros);
  if (init_.isNotNull()) {
    Rcpp::List init(init_);
    // rrr = init["omega"];
    // omega = arma::vec(rrr);
    Rcpp::NumericVector omega_tmp = init["omega"];
    Rcpp::NumericMatrix lambda_tmp = init["lambda"];
    Rcpp::NumericMatrix eta_tmp = init["eta"];
    Rcpp::NumericMatrix phi_tmp = init["phi_mat"];
    Rcpp::NumericMatrix sigmasqeta_tmp = init["prec_eta"];
    Rcpp::NumericMatrix delta_eta_tmp = init["delta_eta"];
    Rcpp::NumericMatrix beta_tmp = init["beta"];
    omega = arma::vec(omega_tmp);
    lambda = Rcpp::as<arma::mat>(lambda_tmp);
    eta = Rcpp::as<arma::mat>(eta_tmp);
    delta_eta = Rcpp::as<arma::mat>(delta_eta_tmp);
    beta = Rcpp::as<arma::mat>(beta_tmp);
    sigmasqeta = Rcpp::as<arma::mat>(sigmasqeta_tmp);
    sigmasqetai = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
    for (arma::uword i = 0; i < dat.nsub; i++) {
      sigmasqetai.rows(i * dat.nreg, (i + 1) * dat.nreg - 1) = sigmasqeta;
    }
    for (arma::uword l = 0; l < dat.ldim; l++) {
      Rcpp::NumericMatrix phi_tmp_tmp =
        phi_tmp(Rcpp::Range(l * dat.nreg, (l + 1) * dat.nreg - 1),
                Rcpp::Range(0, dat.nreg - 1));
      phi.slice(l) = Rcpp::as<arma::mat>(phi_tmp_tmp);
    }
    // rho = init["rho"];
    rho = .1;
    // alpha = init["alpha"];
    alpha = .1;
  } else {
    omega = arma::vec(dat.nreg, arma::fill::ones);
    lambda = arma::mat(dat.basisdim, dat.ldim, arma::fill::randn);
    eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::randn);
    phi = arma::cube(dat.nreg, dat.nreg, dat.ldim, arma::fill::randn);
    sigmasqetai = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
    delta_eta = arma::mat(dat.nreg, dat.ldim, arma::fill::ones);
    beta = arma::mat(dat.designdim * dat.nreg, dat.ldim, arma::fill::randn);
    rho = 0.5;
    alpha = 1;
  }
  compute_sigmasqeta_partial(delta_eta, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
  posterior_omega_shape = dat.nt * dat.nsub / 2. + prior_shape;
  zeta = arma::vec(dat.ldim, arma::fill::ones);
  sigmasqeta = arma::mat(dat.nreg, dat.ldim, arma::fill::ones);
  xi_eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
  delta_beta = arma::mat(dat.nreg * dat.designdim, dat.ldim, arma::fill::ones);
  nu = 60;
  a1 = 2;
  a2 = 2;
  a3 = 2;
}

ParametersWeak::ParametersWeak(const Data& dat, Rcpp::Nullable<Rcpp::List> init_ = R_NilValue) {
  omega_container = arma::mat(dat.nreg, dat.iter);
  zeta_container = arma::mat(dat.ldim, dat.iter);
  lambda_container = arma::cube(dat.basisdim, dat.ldim, dat.iter);
  phi_container = arma::cube(dat.nreg, dat.nreg, dat.iter);
  eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqetai_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqeta_container = arma::cube(dat.nreg, dat.ldim, dat.iter);
  xi_eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  beta_container = arma::cube(dat.designdim * dat.nreg, dat.ldim, dat.iter);
  delta_beta_container = arma::cube(dat.nreg * dat.designdim, dat.ldim, dat.iter);
  delta_eta1_container = arma::cube(dat.nreg, dat.cdim, dat.iter);
  delta_eta2_container = arma::cube(dat.cdim, dat.ldim, dat.iter);
  delta_eta3_container = arma::mat(dat.cdim, dat.iter);
  a1_container = arma::vec(dat.iter);
  a2_container = arma::vec(dat.iter);
  a3_container = arma::vec(dat.iter);
  a4_container = arma::vec(dat.iter);
  nu_container = arma::vec(dat.iter);
  phi = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  delta_beta = arma::mat(dat.nreg * dat.designdim, dat.ldim, arma::fill::ones);
  xi_eta = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
  zeta = arma::vec(dat.ldim, arma::fill::ones);
  sigmasqetai = arma::mat(dat.nsub * dat.nreg, dat.ldim, arma::fill::ones);
  a1 = 2;
  a2 = 2;
  a3 = 2;
  a4 = 2;
  if (init_.isNotNull()) {
    Rcpp::List init(init_);
    Rcpp::NumericMatrix phi_ = init["phi"];
    Rcpp::NumericMatrix beta_ = init["beta"];
    Rcpp::NumericVector omega_ = init["omega"];
    Rcpp::NumericMatrix lambda_ = init["lambda"];
    Rcpp::NumericMatrix delta_eta1_ = init["delta_eta1"];
    Rcpp::NumericMatrix delta_eta2_ = init["delta_eta2"];
    Rcpp::NumericMatrix sigmasqeta_ = init["prec_eta"];
    Rcpp::NumericMatrix eta_ = init["eta"];
    phi = Rcpp::as<arma::mat>(phi_);
    beta = Rcpp::as<arma::mat>(beta_);
    omega = Rcpp::as<arma::vec>(omega_);
    lambda = Rcpp::as<arma::mat>(lambda_);
    delta_eta1 = Rcpp::as<arma::mat>(delta_eta1_);
    delta_eta2 = Rcpp::as<arma::mat>(delta_eta2_);
    // arma::vec tmp1;
    // arma::rowvec tmp2;
    // for (arma::uword c = 0; c < dat.cdim; c++) {
      // tmp1 = delta_eta1.col(c);
      // tmp2 = delta_eta2.row(c);
      // delta_eta1.col(c) = identify_delta1(tmp1);
      // delta_eta2.row(c) = identify_delta2(tmp2);
    // }
    delta_eta1.fill(1);
    delta_eta2.fill(1);
    delta_eta3 = arma::ones(dat.cdim);
    sigmasqeta = Rcpp::as<arma::mat>(sigmasqeta_);
    for (arma::uword i = 0; i < dat.nsub; i++) {
      sigmasqetai.rows(i * dat.nreg, (i + 1) * dat.nreg - 1) = sigmasqeta;
    }
    eta = Rcpp::as<arma::mat>(eta_);
    
  }
  // compute_sigmasqeta_weak(delta_eta1, delta_eta2, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
  nu = 60;
}

void ParametersWeak::update_eta(
    const Data& dat, TransformationsWeak& transf) {
  arma::mat diagomega = arma::diagmat(omega);
  arma::mat diagsigma, Q;
  arma::mat Qbase = arma::kron(arma::eye(dat.ldim, dat.ldim), phi.t() * diagomega * phi);
  arma::mat beta_mat = arma::trans(arma::reshape(beta, dat.designdim, dat.nreg * dat.ldim));
  arma::uword first_eta, last_eta, first, last;
  arma::vec b;
  for (arma::uword i = 0; i < dat.nsub; i++) {
    first_eta = i * dat.nreg;
    last_eta = (i + 1) * dat.nreg - 1;
    first = i * dat.nt;
    last = (i + 1) * dat.nt - 1;
    diagsigma = arma::diagmat(arma::vectorise(sigmasqetai.rows(first_eta, last_eta)));
    if (i == 0) {
      // Rcpp::Rcout << diagsigma << "\n";
    }
    Q = Qbase + diagsigma;
    b = arma::vectorise(phi.t() * diagomega * dat.response.rows(first, last).t() * transf.psi) +
      diagsigma * (beta_mat * dat.design.row(i).t());
    eta.rows(first_eta, last_eta) =
      arma::reshape(bayesreg(b, Q), dat.nreg, dat.ldim);
  }
}

void ParametersWeak::update_lambda(
    const Data& dat, TransformationsWeak& transf) {
  // double eta_sum = 0;
  arma::vec eta_temp;
  arma::mat eta_phi = arma::mat(dat.nsub, dat.nreg);
  arma::mat eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);

  double psi_norm;
  arma::vec eta_vec;
  arma::vec b = arma::zeros(dat.basisdim);
  arma::mat Q, diagomega;
  diagomega = arma::diagmat(omega);
  arma::vec tmprm = arma::vec(dat.basisdim);
  arma::vec tmpad = arma::vec(dat.basisdim);
  for (arma::uword l = 0; l < dat.ldim; l++) {
    // arma::uvec rm = arma::regspace<arma::uvec>(0, dat.ldim - 1);
    // rm.shed_row(l);
    // b.zeros();
    // eta_sum = 0;
    // for (arma::uword i = 0; i < dat.nsub; i++) {
    //   eta_vec = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
    //   eta_sum = eta_sum + arma::dot(eta_vec, eta_vec);
    //   b = b + transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) * diagomega * phi * eta_vec -
    //     transf.btb * lambda.cols(rm) * arma::mat(eta.cols(rm)).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() * phi.t() * diagomega * phi * eta_vec;
    // }
    tmprm.zeros();
    tmpad.zeros();
    b.zeros();
    eta_sum.zeros();
    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_sum = eta_sum + eta_temp * eta_temp.t();
      eta_phi.row(i) = eta_temp.t() * phi.t();
      for (arma::uword lp = 0; lp < dat.ldim; lp++) {
        if (lp != l) {
          tmprm =
            tmprm +
            transf.btb * lambda.col(lp) *
            (eta.col(lp).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() *
            phi.t() * diagomega * eta_phi.row(i).t());
        }
      }
      tmpad =
        tmpad +
        transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) *
        diagomega * eta_phi.row(i).t();
    }
    b = tmpad - tmprm;
    Q = arma::trace(phi.t() * diagomega * phi * eta_sum) * transf.btb + zeta(l) * dat.penalty;
    if (l == 0) {
      lambda.col(l) = bayesreg(b, Q);
    } else {
      lambda.col(l) = bayesreg_orth(b, Q, transf.psi_lin_constr.rows(0, l - 1));
    }
    transf.psi.col(l) = dat.basis * lambda.col(l);
    psi_norm = arma::norm(transf.psi.col(l));
    lambda.col(l) = lambda.col(l) / psi_norm;
    transf.psi.col(l) = transf.psi.col(l) / psi_norm;
    eta.col(l) = eta.col(l) * psi_norm;
    transf.psi_lin_constr.row(l) = transf.psi.col(l).t() * dat.basis;
  }

}

void ParametersWeak::update_omega(
    const Data& dat, TransformationsWeak& transf) {
  double rate, shape;
  shape = .5 * dat.nt * dat.nsub + prior_omega_shape;
  transf.fit.zeros();
  for (arma::uword i = 0; i < dat.nsub; i++) {
    transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.psi *
      eta.rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() * phi.t();
  }
  for (arma::uword r = 0; r < dat.nreg; r++) {
    rate = prior_omega_rate +
      .5 * arma::accu(arma::square(dat.response.col(r) -
    transf.fit.col(r)));
    omega(r) = R::rgamma(shape, 1. / rate);
  }
}

void ParametersWeak::update_phi(
    const Data& dat, TransformationsWeak& transf) {
  double eta_sum = 0;
  double norm;
  arma::vec b = arma::zeros(dat.nreg);
  arma::mat Q = arma::zeros(dat.nreg, dat.nreg);
  arma::mat diagomega = arma::diagmat(omega);
  arma::uword first, last, first_eta, last_eta;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    arma::uvec rm = arma::regspace<arma::uvec>(0, dat.nreg - 1);
    rm.shed_row(r);
    b.zeros();
    eta_sum = 0;
    for (arma::uword i = 0; i < dat.nsub; i++) {
      first = i * dat.nt;
      last = (i + 1) * dat.nt - 1;
      first_eta = i * dat.nreg;
      last_eta = (i + 1) * dat.nreg - 1;
      arma::uvec rm_eta = arma::regspace<arma::uvec>(first_eta, last_eta);
      rm_eta.shed_row(r);
      b = b + arma::vectorise(diagomega * (dat.response.rows(first, last).t() *
        transf.psi * eta.row(first_eta + r).t() - phi.cols(rm) *
        eta.rows(rm_eta) * eta.row(first_eta + r).t()));
      eta_sum = eta_sum + arma::as_scalar(eta.row(first_eta + r) *
        eta.row(first_eta + r).t());
    }
    Q = eta_sum * diagomega + arma::eye(dat.nreg, dat.nreg);
    if (r == 0) {
      phi.col(r) = bayesreg(b, Q);
    } else {
      phi.col(r) = bayesreg_orth(b, Q, transf.phi_lin_constr.rows(0, r - 1));
    }
    norm = arma::norm(phi.col(r));
    phi.col(r) = phi.col(r) / norm;
    transf.phi_lin_constr.row(r) = phi.col(r).t();
    for (arma::uword i = 0; i < dat.nsub; i++) {
      first_eta = i * dat.nreg;
      eta.row(first_eta + r) = eta.row(first_eta + r) * norm;
    }
  }
}

void ParametersWeak::update_a1234(Data& dat) {
  double a_prop;
  double offset = .5;
  double density_old = 0;
  double density_new = 0;
  double mhr;

  a_prop = R::runif(std::max(0., a1 - offset), a1 + offset);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    density_old = R::dgamma(delta_eta1(0, c), a1, 1, true);
    density_new = R::dgamma(delta_eta1(0, c), a_prop, 1, true);
  }
  density_old = density_old + R::dgamma(a1, 2, 1, true);
  density_new = density_new + R::dgamma(a_prop, 2, 1, true);
  mhr = density_new - density_old -
    R::dunif(a_prop, std::max(0., a1 - offset), a1 + offset, 1) +
    R::dunif(a1, std::max(0., a_prop - offset), a_prop + offset, 1);
  if (R::runif(0, 1) < exp(mhr)) {
    a1 = a_prop;
  }

  density_old = 0;
  density_new = 0;
  a_prop = R::runif(std::max(0., a2 - offset), a2 + offset);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 1; r < dat.nreg; r++) {
      density_old = R::dgamma(delta_eta1(r, c), a2, 1, true);
      density_new = R::dgamma(delta_eta1(r, c), a_prop, 1, true);
    }
  }
  density_old = density_old + R::dgamma(a2, 2, 1, true);
  density_new = density_new + R::dgamma(a_prop, 2, 1, true);
  mhr = density_new - density_old -
    R::dunif(a_prop, std::max(0., a2 - offset), a2 + offset, 1) +
    R::dunif(a1, std::max(0., a_prop - offset), a_prop + offset, 1);
  if (R::runif(0, 1) < exp(mhr)) {
    a2 = a_prop;
  }

  density_old = 0;
  density_new = 0;
  a_prop = R::runif(std::max(0., a3 - offset), a3 + offset);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    density_old = R::dgamma(delta_eta2(c, 0), a3, 1, true);
    density_new = R::dgamma(delta_eta2(c, 0), a_prop, 1, true);
  }
  density_old = density_old + R::dgamma(a3, 2, 1, true);
  density_new = density_new + R::dgamma(a_prop, 2, 1, true);
  mhr = density_new - density_old -
    R::dunif(a_prop, std::max(0., a3 - offset), a3 + offset, 1) +
    R::dunif(a3, std::max(0., a_prop - offset), a_prop + offset, 1);
  if (R::runif(0, 1) < exp(mhr)) {
    a3 = a_prop;
  }

  density_old = 0;
  density_new = 0;
  a_prop = R::runif(std::max(0., a4 - offset), a4 + offset);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword l = 1; l < dat.ldim; l++) {
      density_old = R::dgamma(delta_eta2(c, l), a4, 1, true);
      density_new = R::dgamma(delta_eta2(c, l), a_prop, 1, true);
    }
  }
  density_old = density_old + R::dgamma(a4, 2, 1, true);
  density_new = density_new + R::dgamma(a_prop, 2, 1, true);
  mhr = density_new - density_old -
    R::dunif(a_prop, std::max(0., a4 - offset), a4 + offset, 1) +
    R::dunif(a1, std::max(0., a_prop - offset), a_prop + offset, 1);
  if (R::runif(0, 1) < exp(mhr)) {
    a4 = a_prop;
  }
}

void ParametersWeak::update_delta_eta1(Data& dat,
                                       TransformationsWeak& transf) {
  double proposal, density_old, density_new, prior_old, prior_new, p1, p2;
  double M1, M2;
  double step_sigma = .05;
  // double step_sigma;
  double eps = .001;
  arma::mat d1 = delta_eta1;
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      proposal = -1;
      // step_sigma = .2 * arma::stddev(arma::vec(arma::cube(delta_eta1_container.slices(0, current_iter)).tube(r, c)), 1) + eps;
      // if (r == 0 && c == 0) {
        // Rcpp::Rcout << step_sigma << "\n";
      // }
      while (proposal <= 0) {
        proposal = delta_eta1(r, c) + R::rnorm(0, step_sigma);
      }
      d1(r, c) = proposal;
      density_old = get_delta_eta1_density(delta_eta1, delta_eta2,
                                          eta, beta, xi_eta, dat.design, r, c);
      density_new = get_delta_eta1_density(d1, delta_eta2,
                                          eta, beta, xi_eta, dat.design, r, c);
      if (r == 0) {
        prior_old = R::dgamma(delta_eta1(r, c), a1, 1, true);
        prior_new = R::dgamma(proposal, a1, 1, true);

      } else {
        prior_old = R::dgamma(delta_eta1(r, c), a2, 1, true);
        prior_old = R::dgamma(proposal, a2, 1, true);
      }
      p1 = R::pnorm(delta_eta1(r, c), 0, step_sigma, 1, 1);
      p2 = R::pnorm(proposal, 0, step_sigma, 1, 1);
      M1 = density_old + prior_old - p1;
      M2 = density_new + prior_new - p2;
      if (R::runif(0, 1) < exp(M2 - M1)) {
        delta_eta1(r, c) = proposal;
      }
    }
  }
  compute_sigmasqeta_weak(delta_eta1, delta_eta1, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
}

void ParametersWeak::update_delta_eta2(Data& dat,
                                       TransformationsWeak& transf) {
  double proposal, density_old, density_new, prior_old, prior_new, p1, p2;
  // double step_sigma;
  double M1, M2;
  double step_sigma = .05;
  arma::mat d2 = delta_eta2;
  double eps = .001;
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      proposal = -1;
      // step_sigma = .2 * arma::stddev(arma::vec(arma::cube(delta_eta1_container.slices(0, current_iter)).tube(c, l)), 1) + eps;
      while (proposal <= 0) {
        proposal = delta_eta2(c, l) + R::rnorm(0, step_sigma);
      }
      d2(c, l) = proposal;
      density_old = get_delta_eta2_density(delta_eta1, delta_eta2,
                                          eta, beta, xi_eta, dat.design, c, l);
      density_new = get_delta_eta2_density(delta_eta1, d2,
                                          eta, beta, xi_eta, dat.design, c, l);
      if (l == 0) {
        prior_old = R::dgamma(delta_eta2(c, l), a3, 1, true);
        prior_new = R::dgamma(proposal, a3, 1, true);

      } else {
        prior_old = R::dgamma(delta_eta2(c, l), a4, 1, true);
        prior_old = R::dgamma(proposal, a4, 1, true);
      }
      p1 = R::pnorm(delta_eta2(c, l), 0, step_sigma, 1, 1);
      p2 = R::pnorm5(proposal, 0, step_sigma, 1, 1);
      M1 = density_old + prior_old - p1;
      M2 = density_new + prior_new - p2;
      if (R::runif(0, 1) < exp(M2 - M1)) {
        delta_eta2(c, l) = proposal;
      }
    }
  }
  compute_sigmasqeta_weak(delta_eta1, delta_eta1, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
}


void ParametersPartial::update_omega(
    const Data& dat, TransformationsPartial& transf) {
  double rate, shape;
  shape = .5 * dat.nt * dat.nsub;
  arma::vec eta_temp;
  arma::mat eta_phi = arma::mat(dat.nsub, dat.nreg);
  transf.fit.zeros();
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
        transf.psi.col(l) * eta_phi.row(i);
    }
  }
  for (arma::uword r = 0; r < dat.nreg; r++) {
    rate = prior_omega_rate +
        .5 * arma::accu(arma::square(dat.response.col(r) -
        transf.fit.col(r)));
    omega(r) = R::rgamma(shape, 1. / rate);
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

void ParametersPartial::update_eta(const Data& dat, Transformations& transf) {
  arma::uword idx;
  arma::mat eta_sum, eta_phi;
  eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  eta_phi = arma::mat(dat.nsub, dat.nreg);
  arma::uword first, last, first_eta, last_eta;
  arma::vec b, yt, eta_temp;
  arma::mat Q, beta_mat, diagsigma, diagomega = arma::diagmat(omega);
  b = arma::vec(dat.nreg, arma::fill::zeros);
  for (arma::uword i = 0; i < dat.nsub; i++) {
    first = i * dat.nt, last = (i + 1) * dat.nt - 1;
    first_eta = i * dat.nreg, last_eta = (i + 1) * dat.nreg - 1;
    for (arma::uword l = 0; l < dat.ldim; l++) {
      diagsigma = arma::diagmat(sigmasqetai.rows(first_eta, last_eta).col(l));
      beta_mat = arma::reshape(beta.col(l), dat.designdim, dat.nreg).t();
      b = phi.slice(l).t() * diagomega * (dat.response.rows(first, last).t() * transf.psi.col(l)) +
        diagsigma * (beta_mat * dat.design.row(i).t());
      Q = phi.slice(l).t() * diagomega * phi.slice(l) + diagsigma;
      eta.rows(first_eta, last_eta).col(l) = bayesreg(b, Q);
    }
  }
}

void Parameters::update_xi_eta(const Data& dat, Transformations& transf) {
  double shape, rate;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    shape = .5 + nu / 2.;
    arma::uword idx;
    for (arma::uword r = 0; r < dat.nreg; r++) {
      for (arma::uword i = 0; i < dat.nsub; i++) {
        idx = i * dat.nreg + r;
        rate = nu / 2. +
                     .5 * std::pow(eta(idx, l), 2) * sigmasqeta(r, l);
        xi_eta(idx, l) = R::rgamma(shape, 1. / rate);
      }
    }
  }
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
}

void ParametersPartial::update_delta_eta(const Data& dat, Transformations& transf) {
  double tmpsum, ndf, cumprod;
  arma::vec etavec, etavecr, betavec, betavecr, etamean;
  arma::rowvec delta_eta_cumprod_init;

  //////////////

  tmpsum = 0;
  delta_eta(0, 0) = 1;
  compute_sigmasqeta_partial(delta_eta, sigmasqeta);

  for (arma::uword r = 0; r < dat.nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      etavec = eta.col(l);
      etavecr = etavec.rows(r_ind);
      betavec = beta.col(l);
      betavecr = betavec.rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
      etamean = dat.design * betavecr;
      tmpsum = tmpsum + sigmasqeta(r, l) * arma::as_scalar((etavecr - etamean).t() *
        arma::diagmat(arma::mat(xi_eta.col(l)).rows(r_ind)) * (etavecr - etamean));
    }
  }
  ndf = a1 + .5 * dat.nsub * dat.nreg * dat.ldim;
  delta_eta(0, 0) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));

  //////////////

  for (arma::uword l = 1; l < dat.ldim; l++) {
    tmpsum = 0;
    delta_eta(0, l) = 1;
    compute_sigmasqeta_partial(delta_eta, sigmasqeta);

    for (arma::uword rp = 0; rp < dat.nreg; rp++) {
      arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
      for (arma::uword lp = l; lp < dat.ldim; lp++) {
        etavec = eta.col(lp);
        etavecr = etavec.rows(r_ind);
        betavec = beta.col(lp);
        betavecr = betavec.rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
        etamean = dat.design * betavecr;
        tmpsum = tmpsum + sigmasqeta(rp, lp) * arma::as_scalar((etavecr - etamean).t() *
          arma::diagmat(arma::mat(xi_eta.col(lp)).rows(r_ind)) * (etavecr - etamean));
      }
    }
    ndf = a2 + .5 * dat.nsub * (dat.ldim - l) * dat.nreg;
    delta_eta(0, l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
  }

  //////////////

  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 1; r < dat.nreg; r++) {
      tmpsum = 0;
      delta_eta(r, l) = 1;
      compute_sigmasqeta_partial(delta_eta, sigmasqeta);
      for (arma::uword rp = r; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        etavec = eta.col(l);
        etavecr = etavec.rows(r_ind);
        betavec = beta.col(l);
        betavecr = betavec.rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
        etamean = dat.design * betavecr;
        tmpsum = tmpsum + sigmasqeta(rp, l) * arma::as_scalar((etavecr - etamean).t() *
          arma::diagmat(arma::mat(xi_eta.col(l)).rows(r_ind)) * (etavecr - etamean));
      }
      ndf = a3 + .5 * dat.nsub * (dat.nreg - r);
      delta_eta(r, l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
    }
  }
  //////////////

  compute_sigmasqeta_partial(delta_eta, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
}

void ParametersWeak::update_delta_eta_c2(Data& dat) {
  double mult = .05;
  double p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0;
  double proposal;
  arma::mat d1 = delta_eta1;
  arma::mat d2 = delta_eta2;
  arma::vec d3 = delta_eta3;
  double new_sd1, old_sd1, new_sd2, old_sd2;
  double sd3;
  double cut, old_cut;
  double new_density, old_density;
  double mhr;
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      if (r == 0) {
        cut = 1e-6;
        old_cut = 1e-6;
      } else {
        cut = d1(r - 1, c);
        old_cut = delta_eta1(r - 1, c);
      }
      new_sd1 = mult * delta_eta1(r, c);
      proposal = normal_trunc_left(d1(r, c), new_sd1, cut);
      old_sd1 = mult * proposal;
      d1(r, c) = proposal;
      p1 = p1 + R::pnorm5(d1(r, c) - cut, 0, new_sd1, true, true);
      p3 = p3 + R::pnorm5(delta_eta1(r, c) - old_cut, 0, old_sd1, true, true);
    }
    for (arma::uword l = 0; l < dat.ldim; l++) {
      if (l == 0) {
        cut = 1e-6;
        old_cut = 1e-6;
      } else {
        cut = d2(c, l - 1);
        old_cut = delta_eta2(c, l - 1);
      }
      new_sd2 = mult * delta_eta2(c, l);
      proposal = normal_trunc_left(d2(c, l), new_sd2, cut);
      old_sd2 = mult * proposal;
      d2(c, l) = proposal;
      p2 = p1 + R::pnorm5(d2(c, l) - cut, 0, new_sd2, true, true);
      p4 = p3 + R::pnorm5(delta_eta2(c, l) - old_cut, 0, old_sd2, true, true);
    }
    // Rcpp::Rcout << "going in density\n";
    new_density = compute_delta_eta_density_c2(d1, d2, d3, eta, beta,
                                              xi_eta, dat.design, c);
    old_density = compute_delta_eta_density_c2(delta_eta1, delta_eta2, delta_eta3, 
                                               eta, beta, xi_eta, dat.design, c);
    // if (current_iter < .5 * dat.burnin) {
    //   double b = double(current_iter) / double(dat.burnin);
    //   new_density = b * new_density;
    //   old_density = b * old_density;
    // }
    // 
    mhr = new_density - old_density + p3 + p4 - p1 - p2 + p6 - p5;
    // if (c == 0) {
      // Rcpp::Rcout << old_density << "\n" << new_density << "\n" << exp(mhr) << "\n" << "Bob's red mill\n";
    // }
    if (R::runif(0, 1) < exp(mhr)) {
      old_density = new_density;
      delta_eta1.col(c) = d1.col(c);
      delta_eta2.row(c) = d2.row(c);
      delta_eta3(c) = d3(c);
    }
  }
  // Rcpp::Rcout << old_density << "\n";
  // compute_sigmasqeta_weak2(delta_eta1, delta_eta2, sigmasqeta);
  sigmasqeta = delta_eta1 * arma::diagmat(delta_eta3) * delta_eta2;
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
}

void ParametersWeak::update_delta3(Data& dat) {
  arma::vec new_sdv, old_sdv;
  arma::vec d3 = delta_eta3;
  arma::vec etavec, betavec, etamean;
  double mult = .025;
  double cut, old_cut;
  double old_sd, new_sd;
  double new_density = 0, old_density = 0, mhr;
  double p1 = 0, p2 = 0;
  for (arma::uword c = 1; c < dat.cdim; c++) {
    cut = d3(c - 1);
    new_sd = mult * delta_eta3(c);
    d3(c) = normal_trunc_left(d3(c), new_sd, cut);
    old_cut = d3(c - 1);
    old_sd = mult * d3(c);
    new_density = new_density + R::dgamma(d3(c), 1, 1, true) - 
      R::pgamma(cut, 1, 1, false, true);
    old_density = old_density + R::dgamma(delta_eta3(c), 1, 1, true) - 
      R::pgamma(old_cut, 1, 1, false, true);
    p1 = p1 + R::pnorm5(d3(c) - cut, 0, new_sd, true, true);
    p2 = p2 + R::pnorm5(delta_eta3(c) - old_cut, 0, old_sd, true, true);
  }
  arma::mat new_sigmasqeta = delta_eta1 * arma::diagmat(d3) * delta_eta2;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
      etamean = dat.design * betavec;
      old_sdv = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * sigmasqeta(r, l), -.5);
      new_sdv = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * new_sigmasqeta(r, l), -.5);
      old_density = old_density + arma::accu(arma::log_normpdf(etavec, etamean, old_sdv));
      new_density = new_density + arma::accu(arma::log_normpdf(etavec, etamean, new_sdv));
    }
  }

  mhr = new_density - old_density - p2 + p1;
  if (R::runif(0, 1) < exp(mhr)) {
    delta_eta3 = d3;
    sigmasqeta = delta_eta1 * arma::diagmat(delta_eta3) * delta_eta2;
    compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
  }
}

void ParametersWeak::update_delta_eta_c(Data& dat, TransformationsWeak& transf) {
  arma::mat d1 = delta_eta1;
  arma::mat d2 = delta_eta2;
  arma::vec tmp1, sd1;
  arma::rowvec tmp2, sd2;
  double mult = 0.1;
  double old_density, new_density, log_mhr;
  double p1, p2, p3, p4;
  double proposal;
  for (arma::uword c = 0; c < dat.cdim; c++) {
    tmp1 = d1.col(c);
    tmp2 = d2.row(c);
    d1.col(c) = identify_delta1(tmp1);
    d2.row(c) = identify_delta2(tmp2);
    sd1 = mult * d1.col(c);
    sd2 = mult * d2.row(c);
    for (arma::uword r = 0; r < dat.nreg; r++) {
      proposal = d1(r, c) + R::rnorm(0, sd1(r));
      while(proposal < 0) {
        proposal = d1(r, c) + R::rnorm(0, sd1(r));
      }
      d1(r, c) = proposal;
    }
    for (arma::uword l = 0; l < dat.ldim; l++) {
      proposal = d2(c, l) + R::rnorm(0, sd2(l));
      while(proposal + d2(c, l) < 0) {
        proposal = d2(c, l) + R::rnorm(0, sd2(l));
      }
      d2(c, l) = proposal;
    }
    Rcpp::Rcout << d1.col(c) << "\n" << d2.row(c) << "\n";
    new_density = compute_delta_eta_density_c(d1, delta_eta2, eta, beta,
                                xi_eta, dat.design, c, a1, a2, a3, a4);
    old_density = compute_delta_eta_density_c(delta_eta1, delta_eta2, eta, beta,
                                xi_eta, dat.design, c, a1, a2, a3, a4);
    // Rcpp::Rcout << "old_density: " << old_density << "\n" << "new_density: " << new_density << "\n";
    // if (c == 0) {
      // Rcpp::Rcout << d1 << "\n" << delta_eta1 << "\n" << d2 << delta_eta2;
    // }
    p1 = arma::accu(arma::log(arma::normcdf(d1.col(c), arma::zeros(dat.nreg), sd1)));
    p2 = arma::accu(arma::log(arma::normcdf(d2.row(c), arma::zeros<arma::rowvec>(dat.ldim), sd2)));
    p3 = arma::accu(arma::log(arma::normcdf(delta_eta1.col(c), arma::zeros(dat.nreg), sd1)));
    p4 = arma::accu(arma::log(arma::normcdf(delta_eta2.row(c), arma::zeros<arma::rowvec>(dat.ldim), sd2)));
    log_mhr = new_density - old_density + p3 + p4 - p1 - p2;
    // Rcpp::Rcout << "old_density: " << old_density << "\n";
    // Rcpp::Rcout << "new_density: " << new_density << "\n";
    // Rcpp::Rcout << exp(log_mhr) << "\n";
    if (c == 0) {
      // Rcpp::Rcout << old_density << "\n";
    }
    if (R::runif(0, 1) < exp(log_mhr)) {
      // Rcpp::Rcout << "accepted\n";
      delta_eta1.col(c) = d1.col(c);
      delta_eta2.row(c) = d2.row(c);
    }
  }
  compute_sigmasqeta_weak(delta_eta1, delta_eta2, sigmasqeta);
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
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

void ParametersWeak::update_sigmasqeta(const Data& dat) {
  double shape = 1;
  double rate = 1;
  double shape_update = shape + .5 * dat.nsub;
  double rate_update;
  double cut;
  arma::vec etavec, betavec, etamean, diff;
  arma::uvec r_ind;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
      etamean = dat.design * betavec;
      diff = etavec - etamean;
      rate_update = rate + .5 * arma::as_scalar(diff.t() * arma::diagmat(arma::vec(xi_eta.col(l)).rows(r_ind)) * diff);
      if (r == 0 && l == 0) {
        cut = 0;
      } else if (r == 0) {
        cut = sigmasqeta(r, l - 1);
      } else if (l == 0) {
        cut = sigmasqeta(r - 1, l);
      } else {
        cut = std::max(sigmasqeta(r - 1, l), sigmasqeta(r, l - 1));
      }
      sigmasqeta(r, l) = gam_trunc_left(shape_update, 1 / rate_update, cut);
    }
  }
  compute_sigmasqetai(sigmasqeta, xi_eta, sigmasqetai);
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

void ParametersPartial::update_lambda(const Data& dat, Transformations& transf) {
  double psi_norm;
  arma::mat Q, eta_sum, eta_phi, diagomega;
  arma::vec b, eta_temp;
  eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  eta_phi = arma::mat(dat.nsub, dat.nreg);
  diagomega = arma::diagmat(omega);
  b = arma::vec(dat.basisdim);
  arma::mat tmprm = arma::vec(dat.basisdim);
  arma::mat tmpad = arma::vec(dat.basisdim);
  for (arma::uword l = 0; l < dat.ldim; l++) {
    tmprm.zeros();
    tmpad.zeros();
    b.zeros();
    eta_sum.zeros();
    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_sum = eta_sum + eta_temp * eta_temp.t();
      eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      for (arma::uword lp = 0; lp < dat.ldim; lp++) {
        if (lp != l) {
          tmprm =
            tmprm +
            transf.btb * lambda.col(lp) *
            (eta.col(lp).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() *
            phi.slice(lp).t() * diagomega * eta_phi.row(i).t());
        }
      }
      tmpad =
        tmpad +
        transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) *
        diagomega * eta_phi.row(i).t();
    }
    b = tmpad - tmprm;
    Q = arma::trace(phi.slice(l).t() * diagomega * phi.slice(l) * eta_sum) * transf.btb + zeta(l) * dat.penalty;
    if (l == 0) {
      lambda.col(l) = bayesreg(b, Q);
    } else {
    lambda.col(l) = bayesreg_orth(b, Q, transf.psi_lin_constr.rows(0, l - 1));
    }
    transf.psi.col(l) = dat.basis * lambda.col(l);
    psi_norm = arma::norm(transf.psi.col(l));
    lambda.col(l) = lambda.col(l) / psi_norm;
    transf.psi.col(l) = transf.psi.col(l) / psi_norm;
    eta.col(l) = eta.col(l) * psi_norm;
    transf.psi_lin_constr.row(l) = transf.psi.col(l).t() * dat.basis;
  }
}

void Parameters::update_zeta(const Data& dat, Transformations& transf) {
  double shape = prior_zeta_shape + .5 * (dat.penalty_rank);
  double rate;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    rate = prior_zeta_rate + .5 *
      arma::as_scalar(lambda.col(l).t() * dat.penalty * lambda.col(l));
    zeta(l) = R::rgamma(shape, 1.0 / rate);
    zeta(l) = gam_trunc_right(shape, 1.0 / rate, 1e8);
  }
}

void ParametersPartial::update_phi(const Data& dat, TransformationsPartial& transf) {
  double norm;
  arma::cube phi_previous = phi;
  arma::uword first, last, idx;
  arma::vec b, phi_temp;
  arma::uvec r_ind;
  arma::mat Q, diagomega, diageta, eta_sum, C_inv, diag_r;
  b = arma::zeros(dat.nreg * dat.ldim);
  diagomega = arma::diagmat(omega);
  eta_sum = arma::zeros(dat.ldim, dat.ldim);
  // C_inv = arma::inv_sympd(transf.C_rho);
  diag_r = arma::eye(dat.nreg, dat.nreg);
  arma::uvec constr_indices;
  arma::mat tmprm = arma::vec(dat.ldim * dat.nreg);
  arma::mat tmpad = arma::vec(dat.ldim * dat.nreg);
  for (arma::uword r = 0; r < dat.nreg; r++) {
    tmpad.zeros();
    tmprm.zeros();
    eta_sum.zeros();
    r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nreg - 1) * dat.ldim + r);
    for (arma::uword i = 0; i < dat.nsub; i++) {
      first = i * dat.nt, last = (i + 1) * dat.nt - 1;
      idx = i * dat.nreg + r;
      diageta = arma::diagmat(eta.row(idx));
      tmpad = tmpad + arma::vectorise(diagomega * dat.response.rows(first, last).t() *
        transf.psi * diageta);
      for (arma::uword rp = 0; rp < dat.nreg; rp++) {
        if (rp != r) {
          tmprm = tmprm + arma::vectorise(diagomega *
            arma::mat(phi.col(rp)) *
            arma::as_scalar(eta.row(i * dat.nreg + rp) * eta.row(idx).t()));
        }
      }
      eta_sum = eta_sum + diageta * diageta;
    }
    b = tmpad - tmprm;
    Q = arma::kron(eta_sum, diagomega) +
      arma::eye(dat.nreg * dat.ldim, dat.nreg * dat.ldim);
    if (r > 0) {
      constr_indices =
        arma::sort(arma::join_cols(
          constr_indices,
          arma::regspace<arma::uvec>(
            r - 1, dat.nreg, (dat.ldim - 1) * dat.nreg + r)));
    }
    phi_temp = bayesreg_orth(tmpad, Q, transf.phi_lin_constr.rows(constr_indices));
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
  // Align eigenvectors
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
      if (arma::accu(arma::square(phi.slice(l).col(r) + phi_previous.slice(l).col(r))) <
        arma::accu(arma::square(phi.slice(l).col(r) - phi_previous.slice(l).col(r)))) {
        phi.slice(l).col(r) = -phi.slice(l).col(r);
        arma::vec(beta.col(l)).rows(r * dat.designdim,
                  (r + 1) * dat.designdim - 1) =
                    -arma::vec(beta.col(l)).rows(r * dat.designdim,
                               (r + 1) * dat.designdim - 1);
        arma::vec(eta.col(l)).rows(r_ind) = -arma::vec(eta.col(l)).rows(r_ind);
      }
    }
  }
}

void Parameters::update_nu(const Data& dat, Transformations transf) {
  double offset = 10;
  double logratio,
    loglik_old, loglik_new, nu_oldmh, nu_newmh, nu_proposal;
  loglik_old = 0, loglik_new = 0;
  nu_proposal = R::runif(std::max(2., nu - offset), std::min(nu + offset, 128.));
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword i = 0; i < dat.nsub * dat.nreg; i++) {
      loglik_old = loglik_old + R::dgamma(xi_eta(i, l), nu / 2., 2. / nu, true);
      loglik_new = loglik_new + R::dgamma(xi_eta(i, l), nu_proposal / 2., 2. / nu_proposal, true);
    }
    nu_oldmh = loglik_old - R::dunif(nu, std::max(2., nu_proposal - offset),
                                     std::min(nu_proposal + offset, 128.), true);
    nu_newmh = loglik_new - R::dunif(nu_proposal, std::max(2., nu - offset),
                                     std::min(nu + offset, 128.), true);
    logratio = nu_newmh - nu_oldmh;
    if (R::runif(0, 1) < exp(logratio)) {
      nu = nu_proposal;
    }
  }
}

void ParametersPartial::update_a123(const Data& dat) {
  double offset = .5;
  double prior_old, prior_new, logratio, new_logpost,
  loglik_old, loglik_new, a_oldmh, a_newmh, a_proposal;
  a_proposal = R::runif(std::max(0., a1 - offset), a1 + offset);
  loglik_new = 0;
  loglik_old = 0;
  loglik_new = loglik_new + R::dgamma(delta_eta(0, 0), a_proposal, 1, true);
  loglik_old = loglik_old + R::dgamma(delta_eta(0, 0), a1, 1, true);

  prior_new = R::dgamma(a_proposal, 2, 1, true);
  prior_old = R::dgamma(a1, 2, 1, true);
  new_logpost = loglik_new + prior_new;
  old_logpost = loglik_old + prior_old;
  a_newmh = new_logpost - R::dunif(a_proposal, std::max(0., a1 - offset),
                                   a1 + offset, 1);
  a_oldmh = old_logpost - R::dunif(a1, std::max(0., a_proposal - offset),
                                   a_proposal + offset, 1);
  logratio = a_newmh - a_oldmh;
  if (R::runif(0, 1) < exp(logratio)) {
    a1 = a_proposal;
  }

  ////////
  a_proposal = R::runif(std::max(0., a2 - offset), a2 + offset);
  loglik_new = 0; loglik_old = 0;
  for (arma::uword l = 1; l < dat.ldim; l++) {
    loglik_new = loglik_new + R::dgamma(delta_eta.row(0)(l), a_proposal, 1, true);
    loglik_old = loglik_old + R::dgamma(delta_eta.row(0)(l), a2, 1, true);
  }
  prior_new = R::dgamma(a_proposal, 2, 1, true);
  prior_old = R::dgamma(a2, 2, 1, true);
  new_logpost = loglik_new + prior_new;
  old_logpost = loglik_old + prior_old;
  a_newmh = new_logpost - R::dunif(a_proposal, std::max(0., a2 - offset),
                                   a2 + offset, 1);
  a_oldmh = old_logpost - R::dunif(a2, std::max(0., a_proposal - offset),
                                   a_proposal + offset, 1);
  logratio = a_newmh - a_oldmh;
  if (R::runif(0, 1) < exp(logratio)) {
    a2 = a_proposal;
  }

  ////////
  loglik_new = 0; loglik_old = 0;
  a_proposal = R::runif(std::max(0., a3 - offset), a3 + offset);
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 1; r < dat.nreg; r++) {
      loglik_new = loglik_new + R::dgamma(delta_eta(r, l), a_proposal, 1, true);
      loglik_old = loglik_old + R::dgamma(delta_eta(r, l), a3, 1, true);
    }
  }
  prior_new = R::dgamma(a_proposal, 2, 1, true);
  prior_old = R::dgamma(a3, 2, 1, true);
  new_logpost = loglik_new + prior_new;
  old_logpost = loglik_old + prior_old;
  a_newmh = new_logpost - R::dunif(a_proposal, std::max(0., a3 - offset),
                                   a3 + offset, true);
  a_oldmh = old_logpost - R::dunif(a3, std::max(0., a_proposal - offset),
                                   a_proposal + offset, true);
  logratio = a_newmh - a_oldmh;
  if (R::runif(0, 1) < exp(logratio)) {
    a3 = a_proposal;
  }
}

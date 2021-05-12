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
  delta_eta1_container = arma::mat(dat.nreg, dat.iter);
  delta_eta2_container = arma::mat(dat.ldim, dat.iter);
  a1_container = arma::vec(dat.iter);
  a2_container = arma::vec(dat.iter);
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
    Rcpp::NumericVector delta_eta1_ = init["delta_eta1"];
    Rcpp::NumericVector delta_eta2_ = init["delta_eta2"];
    Rcpp::NumericMatrix sigmasqeta_ = init["prec_eta"];
    Rcpp::NumericMatrix eta_ = init["eta"];
    phi = Rcpp::as<arma::mat>(phi_);
    beta = Rcpp::as<arma::mat>(beta_);
    omega = Rcpp::as<arma::vec>(omega_);
    lambda = Rcpp::as<arma::mat>(lambda_);
    delta_eta1 = Rcpp::as<arma::vec>(delta_eta1_);
    delta_eta2 = Rcpp::as<arma::vec>(delta_eta2_);
    sigmasqeta = Rcpp::as<arma::mat>(sigmasqeta_);
    for (arma::uword i = 0; i < dat.nsub; i++) {
      sigmasqetai.rows(i * dat.nreg, (i + 1) * dat.nreg - 1) = sigmasqeta;
    }
    eta = Rcpp::as<arma::mat>(eta_);
  }
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
    Q = Qbase + diagsigma;
    b = arma::vectorise(phi.t() * diagomega * dat.response.rows(first, last).t() * transf.psi) +
      diagsigma * (beta_mat * dat.design.row(i).t());
    eta.rows(first_eta, last_eta) = 
      arma::reshape(bayesreg(b, Q), dat.nreg, dat.ldim);
  }
}

void ParametersWeak::update_lambda(
    const Data& dat, TransformationsWeak& transf) {
  double eta_sum = 0;
  double psi_norm;
  arma::vec eta_vec;
  arma::vec b = arma::zeros(dat.basisdim);
  arma::mat Q, diagomega;
  diagomega = arma::diagmat(omega);
  for (arma::uword l = 0; l < dat.ldim; l++) {
    arma::uvec rm = arma::regspace<arma::uvec>(0, dat.ldim - 1);
    rm.shed_row(l);
    b.zeros();
    eta_sum = 0;
    for (arma::uword i = 0; i < dat.nsub; i++) {
      eta_vec = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_sum = eta_sum + arma::dot(eta_vec, eta_vec);
      b = b + transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) * diagomega * phi * eta_vec - 
        transf.btb * lambda.cols(rm) * arma::mat(eta.cols(rm)).rows(i * dat.nreg, (i + 1) * dat.nreg - 1).t() * phi.t() * diagomega * phi * eta_vec;
    }
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

void ParametersWeak::update_a123(const Data& dat) {
  
}

void ParametersWeak::update_delta_eta1(const Data& dat,
                                       TransformationsWeak& transf) {
  double ndf, tmpsum;
  arma::vec etavec, betavec, etamean;
  arma::mat delta_eta_cumprod;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    tmpsum = 0;
    delta_eta1(r) = 1;
    delta_eta_cumprod = arma::cumprod(delta_eta1) * 
      arma::cumprod(delta_eta2).t();
    for (arma::uword rp = r; rp < dat.nreg; rp++) {
      arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
      for (arma::uword l = 0; l < dat.ldim; l++) {
        etavec = arma::vec(eta.col(l)).rows(r_ind);
        betavec = arma::vec(beta.col(l)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
        etamean = dat.design * betavec;
        tmpsum = tmpsum + delta_eta_cumprod(rp, l) * arma::as_scalar((etavec - etamean).t() *
          arma::diagmat(arma::mat(xi_eta.col(l)).rows(r_ind)) * (etavec - etamean));
      }
    }
    if (r == 0) {
      ndf = a1 + .5 * dat.nsub * dat.nreg * dat.ldim;
      delta_eta1(r) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));

    } else {
      ndf = a2 + .5 * dat.nsub * (dat.nreg - r) * dat.ldim;
      delta_eta1(r) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
    }
  }
  transf.delta_eta_cumprod = arma::cumprod(delta_eta1) * 
    arma::cumprod(delta_eta2).t();
}

void ParametersWeak::update_delta_eta2(const Data& dat,
                                       TransformationsWeak& transf) {
  double ndf, tmpsum;
  arma::vec etavec, betavec, etamean;
  arma::mat delta_eta_cumprod;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    tmpsum = 0;
    delta_eta2(l) = 1;
    delta_eta_cumprod = arma::cumprod(delta_eta1) * 
      arma::cumprod(delta_eta2).t();
    for (arma::uword r = 0; r < dat.nreg; r++) {
      arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
      for (arma::uword lp = l; lp < dat.ldim; lp++) {
        etavec = arma::vec(eta.col(lp)).rows(r_ind);
        betavec = arma::vec(beta.col(lp)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
        etamean = dat.design * betavec;
        tmpsum = tmpsum + delta_eta_cumprod(r, lp) * arma::as_scalar((etavec - etamean).t() *
          arma::diagmat(arma::mat(xi_eta.col(lp)).rows(r_ind)) * (etavec - etamean));
      }
    }
    if (l == 0) {
      ndf = a3 + .5 * dat.nsub * dat.nreg * dat.ldim;
      delta_eta2(l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
    } else {
      ndf = a4 + .5 * dat.nsub * dat.nreg * (dat.ldim - l);
      delta_eta2(l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
    }
  }
  transf.delta_eta_cumprod = arma::cumprod(delta_eta1) * 
    arma::cumprod(delta_eta2).t();
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
  // transf.fit.zeros();
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      for (arma::uword r = 0; r < dat.nreg; r++) {
        idx = i * dat.nreg + r;
        sigmasqetai(idx, l) = xi_eta(idx, l) * transf.delta_eta_cumprod(r, l);
      }
    }
  }
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
                     .5 * std::pow(eta(idx, l), 2) * transf.delta_eta_cumprod(r, l);
        xi_eta(idx, l) = R::rgamma(shape, 1. / rate);
      }
    }
  }
}

void ParametersPartial::update_delta_eta(const Data& dat, Transformations& transf) {
  double tmpsum, ndf;
  arma::mat delta_eta_cumprod;
  arma::vec etavec, etavecr, betavec, betavecr, etamean;
  delta_eta_cumprod = arma::mat(dat.nreg, dat.ldim);
  arma::rowvec delta_eta_cumprod_init;
  
  //////////////

  tmpsum = 0;
  delta_eta(0, 0) = 1;
  delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
  for (arma::uword l = 0; l < dat.ldim; l++) {
    delta_eta_cumprod.col(l) = arma::cumprod(delta_eta.col(l));
    if (l > 0) {
      delta_eta_cumprod.col(l) = delta_eta_cumprod.col(l) *
        delta_eta_cumprod_init(l - 1);
    }
  }
  for (arma::uword r = 0; r < dat.nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      etavec = eta.col(l);
      etavecr = etavec.rows(r_ind);
      betavec = beta.col(l);
      betavecr = betavec.rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
      etamean = dat.design * betavecr;
      tmpsum = tmpsum + delta_eta_cumprod(r, l) * arma::as_scalar((etavecr - etamean).t() *
        arma::diagmat(arma::mat(xi_eta.col(l)).rows(r_ind)) * (etavecr - etamean));
    }
  }
  ndf = a1 + .5 * dat.nsub * dat.nreg * dat.ldim;
  delta_eta(0, 0) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
  
  //////////////
  
  for (arma::uword l = 1; l < dat.ldim; l++) {
    tmpsum = 0;
    delta_eta(0, l) = 1;
  delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
    for (arma::uword lp = 0; lp < dat.ldim; lp++) {
      delta_eta_cumprod.col(lp) = arma::cumprod(delta_eta.col(lp));
      if (lp > 0) {
        delta_eta_cumprod.col(lp) = delta_eta_cumprod.col(lp) *
          delta_eta_cumprod_init(lp - 1);
      }
    }
    for (arma::uword rp = 0; rp < dat.nreg; rp++) {
      arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
      for (arma::uword lp = l; lp < dat.ldim; lp++) {
        etavec = eta.col(lp);
        etavecr = etavec.rows(r_ind);
        betavec = beta.col(lp);
        betavecr = betavec.rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
        etamean = dat.design * betavecr;
        tmpsum = tmpsum + delta_eta_cumprod(rp, lp) * arma::as_scalar((etavecr - etamean).t() *
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
      delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
      delta_eta_cumprod.col(0) = arma::cumprod(delta_eta.col(0));
      for (arma::uword lp = 1; lp < dat.ldim; lp++) {
        delta_eta_cumprod.col(lp) = arma::cumprod(delta_eta.col(lp)) *
          delta_eta_cumprod_init(lp - 1);
      }
      for (arma::uword rp = r; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        etavec = eta.col(l);
        etavecr = etavec.rows(r_ind);
        betavec = beta.col(l);
        betavecr = betavec.rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
        etamean = dat.design * betavecr;
        tmpsum = tmpsum + delta_eta_cumprod(rp, l) * arma::as_scalar((etavecr - etamean).t() *
          arma::diagmat(arma::mat(xi_eta.col(l)).rows(r_ind)) * (etavecr - etamean));
      }
      ndf = a3 + .5 * dat.nsub * (dat.nreg - r);
      delta_eta(r, l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
    }
  }
  
  //////////////
  
  delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
  transf.delta_eta_cumprod.col(0) = arma::cumprod(delta_eta.col(0));
  for (arma::uword l = 1; l < dat.ldim; l++) {
    transf.delta_eta_cumprod.col(l) = arma::cumprod(delta_eta.col(l)) * 
      delta_eta_cumprod_init(l - 1);
  }
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
  transf.fit.zeros();
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
        arma::vec(transf.fit_eta.col(l)).rows(r_ind) = -arma::vec(transf.fit_eta.col(l)).rows(r_ind);
      }
    }
  }
}

void Parameters::update_nu(const Data& dat, Transformations transf) {
  double offset = 2;
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

#include "parameters.h"

Parameters::Parameters(Data& dat, Rcpp::Nullable<Rcpp::List> init_ = R_NilValue) {
  omega_container = arma::mat(dat.nreg, dat.iter);
  zeta_container = arma::mat(dat.ldim, dat.iter);
  lambda_container = arma::cube(dat.basisdim, dat.ldim, dat.iter);
  phi_container = arma::field<arma::cube>(dat.iter);
  phi0_container = arma::cube(dat.nreg, dat.nreg, dat.iter);
  eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqetai_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  sigmasqeta_container = arma::cube(dat.nreg, dat.ldim, dat.iter);
  xi_eta_container = arma::cube(dat.nsub * dat.nreg, dat.ldim, dat.iter);
  beta_container = arma::cube(dat.designdim * dat.nreg, dat.ldim, dat.iter);
  delta_beta_container = arma::cube(dat.nreg * dat.designdim, dat.ldim, dat.iter);
  delta_eta_container = arma::cube(dat.nreg, dat.ldim, dat.iter);
  tau_phi0_container = arma::cube(dat.nreg, dat.nreg, dat.iter);
  rho_container = arma::vec(dat.iter);
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
    Rcpp::NumericMatrix phi_tmp = init["phi"];
    Rcpp::NumericMatrix sigmasqeta_tmp = init["preceta"];
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
    init_phi = phi;
    rho = init["rho"];
    alpha = init["alpha"];

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
  phi0 = arma::mat(dat.nreg, dat.nreg);
  for (arma::uword r = 0; r < dat.nreg; r++) {
    phi0.col(r) = arma::vec(arma::mean(phi.col(r), 2));
  }
  tau_phi0 = arma::mat(dat.nreg, dat.nreg, arma::fill::ones);
  a1 = 2;
  a2 = 2;
  a3 = 2;
}

void Parameters::update_omega(Data& dat, Transformations& transf) {
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

void Parameters::update_eta(Data& dat, Transformations& transf) {
  arma::mat eta_sum, eta_phi;
  eta_sum = arma::mat(dat.nreg, dat.nreg, arma::fill::zeros);
  eta_phi = arma::mat(dat.nsub, dat.nreg);
  arma::uword first, last, first_eta, last_eta;
  arma::vec b, yt, eta_temp;
  arma::mat Q, beta_mat, diagsigma, diagomega = arma::diagmat(omega);
  // transf.fit.zeros();
  b = arma::vec(dat.nreg, arma::fill::zeros);
  for (arma::uword i = 0; i < dat.nsub; i++) {
    first = i * dat.nt, last = (i + 1) * dat.nt - 1;
    first_eta = i * dat.nreg, last_eta = (i + 1) * dat.nreg - 1;
    for (arma::uword l = 0; l < dat.ldim; l++) {
      diagsigma = arma::diagmat(sigmasqetai.rows(first_eta, last_eta).col(l));
      // yt = arma::trans((transf.psi.col(l).t() * dat.response.rows(first, last)) * phi.slice(l));
      beta_mat = arma::reshape(beta.col(l), dat.designdim, dat.nreg).t();
      // b = phi.slice(l).t() * diagomega * phi.slice(l) * yt +
        // diagsigma * (beta_mat * dat.design.row(i).t());
      b = phi.slice(l).t() * diagomega * (dat.response.rows(first, last).t() * transf.psi.col(l)) +
        diagsigma * (beta_mat * dat.design.row(i).t());
      Q = phi.slice(l).t() * diagomega * phi.slice(l) + diagsigma;
      eta.rows(first_eta, last_eta).col(l) = bayesreg(b, Q);
      // transf.fit.rows(first, last) = transf.fit.rows(first, last) +
        // transf.psi.col(l) * (eta.rows(first_eta, last_eta).col(l).t() * phi.slice(l).t());
    }
  }
  // transf.fit.zeros();
  // for (arma::uword i = 0; i < dat.nsub; i++) {
  //   for (arma::uword l = 0; l < dat.ldim; l++) {
  //     eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
  //     eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
  //     transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
  //       transf.psi.col(l) * eta_phi.row(i);
  //   }
  // }
}

void Parameters::update_xi_eta(Data& dat, Transformations& transf) {
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

void Parameters::update_delta_eta(Data& dat, Transformations& transf) {
  // Rcpp::Rcout << 1 / delta_eta(0,0) << "\n";
  double tmpsum, ndf, cumprod;
  arma::mat delta_eta_cumprod;
  arma::vec etavec, etavecr, betavec, betavecr, etamean;
  delta_eta_cumprod = arma::mat(dat.nreg, dat.ldim);
  arma::rowvec delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
  // Rcpp::Rcout << "sigmasqetai \n" << sigmasqetai.rows(0, 4) << "\n";
  for (arma::uword l = 0; l < dat.ldim; l++) {
    delta_eta_cumprod.col(l) = arma::cumprod(delta_eta.col(l));
    if (l > 0) {
      delta_eta_cumprod.col(l) = delta_eta_cumprod.col(l) * 
        delta_eta_cumprod_init(l - 1);
    }
  }
  arma::uword l = 0;
  arma::uword r = 0;
  // Rcpp::Rcout << "delta_eta\n " << delta_eta << "\n"; 
  // Rcpp::Rcout << "new iter: \n" << delta_eta_cumprod << "\n";
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      tmpsum = 0;
      delta_eta(r, l) = 1;
      delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
      delta_eta_cumprod.col(0) = arma::cumprod(delta_eta.col(0));
      for (arma::uword lp = 1; lp < dat.ldim; lp++) {
        delta_eta_cumprod.col(lp) = arma::cumprod(delta_eta.col(lp)) * 
          delta_eta_cumprod_init(lp - 1);
      }
      
      for (arma::uword lp = 0; lp < l; lp++) {
        delta_eta_cumprod.col(lp).zeros();
      }
      for (arma::uword rp = 0; rp < r; rp++) {
        delta_eta_cumprod.col(l).row(rp) = 0;
      }
      
      if (r > 0) {
        for (arma::uword lp = 0; lp < dat.ldim; lp++) {
          if (lp != l) delta_eta_cumprod.col(lp).zeros();
        }
      }
      // if (r == 0 & l == 0) {
        
        // Rcpp::Rcout << "this is delta_eta\n" << delta_eta << "\n";
        // Rcpp::Rcout << "this is delta_eta_cumprod \n" << delta_eta_cumprod << "\n";
      // }
      for (arma::uword rp = 0; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = 0; lp < dat.ldim; lp++) {
          arma::uvec lpv = arma::ones<arma::uvec>(1) * lp;
          etavec = eta.col(lp);
          etavecr = etavec.rows(r_ind);
          betavec = beta.col(lp);
          betavecr = betavec.rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavecr;
          // xi_eta.submat(r_ind, lpv)
          // tmpsum = tmpsum + delta_eta_cumprod(rp, lp) * 
          //   arma::as_scalar((etavecr - etamean).t() *
          //   arma::diagmat(xi_eta.submat(r_ind, lpv)) * 
          //   (etavecr - etamean));
          tmpsum = tmpsum + delta_eta_cumprod(rp, lp) * arma::as_scalar((etavecr - etamean).t() *
            arma::diagmat(xi_eta.submat(r_ind, lpv)) * (etavecr - etamean));
          if (r == 0) {
            if (l == 0) {
              //     Rcpp::Rcout << etavecr.rows(0, 5) << "\n";
              // Rcpp::Rcout << "beta = " << 0 << "rp = " << rp << "  lp = " << lp << "  mutlply_this = " << delta_eta_cumprod(rp, lp)  << "  add this = " << arma::as_scalar((etavecr - etamean).t() * (etavecr - etamean)) << "\n";
              // Rcpp::Rcout << arma::as_scalar((etavecr).t() * xi_eta.submat(r_ind, lpv)) * (etavecr)) << "\n";
            }
          }
            
            // tmpsum = tmpsum + 
            //   arma::as_scalar((etavecr).t() * 
            //   arma::diagmat(arma::vectorise(delta_eta_cumprod(rp, lp) * xi_eta.submat(r_ind, lpv))) *
            //   (etavecr));
          }
        }
        if (l == 0 & r == 0) {
          ndf = a1 + .5 * dat.nsub * dat.nreg * dat.ldim;
        } else if (r == 0 & l > 0) {
          ndf = a2 + .5 * dat.nsub * dat.nreg * (dat.ldim - l);
        } else {
          ndf = a3 + .5 * dat.nsub * (dat.nreg - r);
        }
        if (r == 0 & l == 0) {
          // Rcpp::Rcout << "tmpsum = " << tmpsum << "\n";
        }
        // Rcpp::Rcout << tmpsum << "\n";
        delta_eta(r, l) = R::rgamma(ndf, 1. / (1 + .5 * tmpsum));
        // Rcpp::Rcout << eta.rows(0, 3) << "\n";
    }
  }
  delta_eta_cumprod_init = arma::cumprod(delta_eta.row(0));
  delta_eta_cumprod.col(0) = arma::cumprod(delta_eta.col(0));
  for (arma::uword l = 1; l < dat.ldim; l++) {
    delta_eta_cumprod.col(l) = arma::cumprod(delta_eta.col(l)) * 
      delta_eta_cumprod_init(l - 1);
  }
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword i = 0; i < dat.nsub; i++) {
      sigmasqetai.submat(dat.nreg * i, l, (i + 1) * dat.nreg - 1, l) = 
        xi_eta.submat(i * dat.nreg, l, (i + 1) * dat.nreg - 1, l) %
        delta_eta_cumprod.col(l);
    }
  }
}

void Parameters::update_delta_eta_mh(const Data& dat, Transformations& transf) {
  double loglik_old, loglik_new, prior_old, prior_new, logpost_new, logpost_old;
  double offset = .5;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      loglik_old = 0;
      loglik_new = 0;
      
    }
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

void Parameters::update_lambda(const Data& dat, Transformations& transf) {
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
      // eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      // eta_sum = eta_sum + eta_temp * eta_temp.t();
      // eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      // transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) -
        // transf.psi.col(l) * eta_phi.row(i);
      // b = tmpad - tmprm;
      // b = b + (transf.bty.rows(i * dat.basisdim, (i + 1) * dat.basisdim - 1) -
        // dat.basis.t() * transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1)) *
        // diagomega * (eta_phi.row(i)).t();
  // }
    Q = arma::trace(phi.slice(l).t() * diagomega * phi.slice(l) * eta_sum) * transf.btb + zeta(l) * dat.penalty;
    // arma::uvec indices =arma::regspace<arma::uvec>(0, 1, dat.ldim - 1);
    // indices.shed_row(l);
    // lambda.col(l) = bayesreg(b, Q);
    if (l == 0) {
      lambda.col(l) = bayesreg(b, Q);
    } else {
    lambda.col(l) = bayesreg_orth(b, Q, transf.psi_lin_constr.rows(0, l - 1));
    }
    // psi_norm = 1;
    transf.psi.col(l) = dat.basis * lambda.col(l);
    psi_norm = arma::norm(transf.psi.col(l));
    lambda.col(l) = lambda.col(l) / psi_norm;
    transf.psi.col(l) = transf.psi.col(l) / psi_norm;
    eta.col(l) = eta.col(l) * psi_norm;
    transf.psi_lin_constr.row(l) = transf.psi.col(l).t() * dat.basis;
    // for (arma::uword i = 0; i < dat.nsub; i++) {
      // eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      // eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      // transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
        // transf.psi.col(l) * eta_phi.row(i);
    // }
  }
  transf.fit.zeros();
  for (arma::uword i = 0; i < dat.nsub; i++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      eta_temp = eta.col(l).rows(i * dat.nreg, (i + 1) * dat.nreg - 1);
      eta_phi.row(i) = eta_temp.t() * phi.slice(l).t();
      transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) = transf.fit.rows(i * dat.nt, (i + 1) * dat.nt - 1) +
        transf.psi.col(l) * eta_phi.row(i);
    }
  }
}

void Parameters::update_zeta(const Data& dat, Transformations& transf) {
  double shape = prior_zeta_shape + .5 * (dat.basisdim - 2);
  double rate;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    rate = prior_zeta_rate + .5 *
      arma::as_scalar(lambda.col(l).t() * dat.penalty * lambda.col(l));
    zeta(l) = R::rgamma(shape, 1.0 / rate);
  }
}

void Parameters::update_phi(const Data& dat, Transformations& transf) {
  double norm;
  arma::uword first, last, idx;
  arma::vec b, phi_temp;
  arma::uvec r_ind;
  arma::mat Q, diagomega, diageta, eta_sum, C_inv, diag_r;
  b = arma::zeros(dat.nreg * dat.ldim);
  diagomega = arma::diagmat(omega);
  eta_sum = arma::zeros(dat.ldim, dat.ldim);
  C_inv = arma::inv_sympd(transf.C_rho);
  diag_r = arma::eye(dat.nreg, dat.nreg);
  transf.fit.zeros();
  arma::uvec constr_indices;
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
    b = b + arma::vectorise(arma::reshape(phi0.col(r), dat.nreg, dat.ldim) * C_inv);
    Q = arma::kron(eta_sum, diagomega) +
      arma::kron(C_inv, diag_r);
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
  // Align eigenvectors
  // for (arma::uword l = 0; l < dat.ldim; l++) {
  //   for (arma::uword r = 0; r < dat.nreg; r++) {
  //     arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
  //     if (arma::accu(arma::square(phi.slice(l).col(r) + init_phi.slice(l).col(r))) <
  //       arma::accu(arma::square(phi.slice(l).col(r) - init_phi.slice(l).col(r)))) {
  //       Rcpp::Rcout << "l = " << l << "  r = " << r << "\n";
  //       phi.slice(l).col(r) = -phi.slice(l).col(r);
  //       // arma::vec(beta.col(l)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1) = -arma::vec(beta.col(l)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
  //       // arma::vec(eta.col(l)).rows(r_ind) = -arma::vec(eta.col(l)).rows(r_ind);
  //       // arma::vec(transf.fit_eta.col(l)).rows(r_ind) = -arma::vec(transf.fit_eta.col(l)).rows(r_ind);
  //     }
  //   }
  // }
}

void Parameters::update_phi0(const Data& dat, Transformations &transf) {
  arma::vec ones_v = arma::ones(dat.ldim);
  arma::mat C_inv;
  double Q, b;
  C_inv = arma::inv_sympd(transf.C_rho);
  Q = arma::as_scalar(
    ones_v.t() * arma::inv_sympd(transf.C_rho) * ones_v);
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword rp = 0; rp < dat.nreg; rp++) {
      
      b = arma::as_scalar(ones_v.t() * C_inv * arma::vec(phi.tube(rp, r)));
      phi0(rp, r) = R::rnorm(b / Q, std::pow(Q, -.5));
      if (rp == 0) {
        if (r == 0) {
          Rcpp::Rcout << "phi tube: " << phi.tube(rp, r) <<"\n";
          Rcpp::Rcout << "mean: " << b / Q <<"\n";
          Rcpp::Rcout << "sd: " << std::pow(Q, -.5) << "\n";
        }
      }
    }
  }
}

void Parameters::update_tau_phi0(const Data& dat, Transformations &transf) {
  
}
void Parameters::update_rho(const Data &dat, Transformations &transf) {
  double offset = .025;
  double prior_old, prior_new, logratio, new_logpost,
    loglik_old, loglik_new, rho_oldmh, rho_newmh, rho_proposal;
  rho_proposal = R::runif(std::max(0., rho - offset), std::min(rho + offset, 1.));
  arma::mat C_rho_proposal = alpha * rho_proposal * transf.ones_mat + 
                             alpha * (1 - rho_proposal) * arma::eye(dat.ldim, dat.ldim);
  arma::mat x = arma::mat(dat.nreg * dat.nreg, dat.ldim);
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword rr = 0; rr < dat.nreg; rr++) {
      x.row(rr + r * dat.nreg) = arma::vec(phi.tube(rr, r)).t();
    }
  }
  loglik_old = arma::accu(dmvnrm_arma_fast(
    x, arma::zeros<arma::rowvec>(dat.ldim), transf.C_rho, true));
  
  loglik_new = arma::accu(dmvnrm_arma_fast(
    x, arma::zeros<arma::rowvec>(dat.ldim), C_rho_proposal, true));
  
  prior_old = R::dbeta(rho, rho_shape1, rho_shape2, true);
  prior_new = R::dbeta(rho_proposal, rho_shape1, rho_shape2, true);
  old_logpost = loglik_old + prior_old;
  new_logpost = loglik_new + prior_new;
  rho_oldmh = old_logpost - R::dunif(rho, std::max(0., rho_proposal - offset),
                                     std::min(rho_proposal + offset, 1.), true);
  rho_newmh = new_logpost - R::dunif(rho_proposal, std::max(0., rho - offset),
                                     std::min(rho + offset, 1.), true);
  logratio = rho_newmh - rho_oldmh;
  if (R::runif(0, 1) < exp(logratio)) {
    rho = rho_proposal;
    transf.C_rho = C_rho_proposal;
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

void Parameters::update_a123(const Data& dat, Transformations& transf) {
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

void Parameters::update_alpha(const Data& dat, Transformations& transf) {
  arma::mat C_rho_inv;
  double shape, rate;
  shape = .5 + .5 * dat.nreg * dat.nreg * dat.ldim;
  rate = .5;
  C_rho_inv = arma::inv_sympd(rho * transf.ones_mat + 
    (1 - rho) * arma::eye(dat.ldim, dat.ldim));
  for (arma::uword r = 0; r < dat.nreg; r++) {
    for (arma::uword rp = 0; rp < dat.nreg; rp++) {
      rate = rate + .5 * arma::as_scalar(arma::vec(phi.tube(rp, r)).t() * 
        C_rho_inv * 
        arma::vec(phi.tube(rp, r)));
    }
  }
  alpha = 1. / R::rgamma(shape, 1. / (1 + rate));
}
void Parameters::update_a12(const Data& dat, Transformations& transf) {
  double offset = .5;
  double prior_old, prior_new, logratio, new_logpost,
  loglik_old, loglik_new, a_oldmh, a_newmh, a_proposal;
  a_proposal = R::runif(std::max(0., a1 - offset), a1 + offset);
  loglik_new = 0;
  loglik_old = 0;
  for (arma::uword l = 0; l < dat.ldim; l++) {
    loglik_new = loglik_new + R::dgamma(delta_eta.row(0)(l), a_proposal, 1, true);
    loglik_old = loglik_old + R::dgamma(delta_eta.row(0)(l), a1, 1, true);
  }
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
  loglik_new = 0; loglik_old = 0;
  a_proposal = R::runif(std::max(0., a2 - offset), a2 + offset);
  for (arma::uword l = 1; l < dat.ldim; l++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      loglik_new = loglik_new + R::dgamma(delta_eta(r, l), a_proposal, 1, true);
      loglik_old = loglik_old + R::dgamma(delta_eta(r, l), a2, 1, true);
    }
  }
  prior_new = R::dgamma(a_proposal, 2, 1, true);
  prior_old = R::dgamma(a2, 2, 1, true);
  new_logpost = loglik_new + prior_new;
  old_logpost = loglik_old + prior_old;
  a_newmh = new_logpost - R::dunif(a_proposal, std::max(0., a2 - offset),
                                   a2 + offset, true);
  a_oldmh = old_logpost - R::dunif(a2, std::max(0., a_proposal - offset), 
                                   a_proposal + offset, true);
  logratio = a_newmh - a_oldmh;
  if (R::runif(0, 1) < exp(logratio)) {
    a2 = a_proposal;
  }
}
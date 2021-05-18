#include "sampler.h"

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List Sampler::write_data() {
  return(dat.write_data());
}

Rcpp::List Sampler::write_control() {
  return(Rcpp::List::create(
      Rcpp::Named("iterations", dat.iter),
      Rcpp::Named("thin", dat.thin),
      Rcpp::Named("burnin", dat.burnin)
  ));
}

SamplerPartial::SamplerPartial(Data& dat_, Rcpp::Nullable<Rcpp::List> init_) {
  dat = dat_;
  pars = ParametersPartial(dat, init_);
  transf = TransformationsPartial(dat, pars);
}

SamplerWeak::SamplerWeak(Data& dat_, Rcpp::Nullable<Rcpp::List> init_) {
  dat = dat_;
  pars = ParametersWeak(dat, init_);
  transf = TransformationsWeak(dat, pars);
  arma::uword dim = pars.delta_eta1.n_elem + pars.delta_eta2.n_elem;
  p = arma::randn(dim);
  q = arma::vec(dim);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      q.row(c * dat.nreg + r) = pars.delta_eta1.col(c).row(r);
    }
  }
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      q.row(dat.cdim * dat.nreg + c * dat.ldim + l) = 
        pars.delta_eta2.row(c).col(l);
    }
  }
  
}

void SamplerPartial::sample() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort()) {
        Rcpp::Rcout << "MCMC terminated by user\n";
        goto stop;
      }
      transf.complete_response(dat, pars);
      pars.update_lambda(dat, transf);
      pars.update_zeta(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_eta(dat, transf);
      pars.update_xi_eta(dat, transf);
      pars.update_delta_eta(dat, transf);
      pars.update_beta(dat, transf);
      pars.update_delta_beta(dat, transf);
      pars.update_omega(dat, transf);
      pars.update_nu(dat, transf);
      pars.update_a123(dat);
    }
    progress_bar.increment();
    write_samples();
  }
  stop:
    NULL;
}

void SamplerPartial::write_samples() {
  pars.lambda_container.slice(pars.current_iter) = pars.lambda;
  pars.beta_container.slice(pars.current_iter) = pars.beta;
  pars.delta_beta_container.slice(pars.current_iter) = pars.delta_beta;
  pars.delta_eta_container.slice(pars.current_iter) = pars.delta_eta;
  pars.omega_container.col(pars.current_iter) = pars.omega;
  pars.xi_eta_container.slice(pars.current_iter) = pars.xi_eta;
  pars.zeta_container.col(pars.current_iter) = pars.zeta;
  pars.eta_container.slice(pars.current_iter) = pars.eta;
  pars.phi_container(pars.current_iter) = pars.phi;
  pars.sigmasqetai_container.slice(pars.current_iter) = pars.sigmasqetai;
  pars.delta_eta_container.slice(pars.current_iter) = pars.delta_eta;
  pars.nu_container(pars.current_iter) = pars.nu;
  pars.a1_container(pars.current_iter) = pars.a1;
  pars.a2_container(pars.current_iter) = pars.a2;
  pars.a3_container(pars.current_iter) = pars.a3;
  current_iter++;
}

Rcpp::List SamplerPartial::get_samples() {
  return Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("delta_beta", pars.delta_beta_container),
                            Rcpp::Named("delta_eta", pars.delta_eta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("omega", pars.omega_container),
                            Rcpp::Named("xi_eta", pars.xi_eta_container),
                            Rcpp::Named("zeta", pars.zeta_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("sigmasqetai", pars.sigmasqetai_container),
                            Rcpp::Named("delta_eta", pars.delta_eta_container),
                            Rcpp::Named("nu", pars.nu_container),
                            Rcpp::Named("a1", pars.a1_container),
                            Rcpp::Named("a2", pars.a2_container),
                            Rcpp::Named("a3", pars.a3_container),
                            Rcpp::Named("fit", transf.fit));
}

void SamplerWeak::sample() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort()) {
        Rcpp::Rcout << "MCMC terminated by user\n";
        goto stop;
      }
      pars.beta.zeros();
      transf.complete_response(dat, pars);
      pars.update_lambda(dat, transf);
      pars.update_zeta(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_eta(dat, transf);
      pars.update_xi_eta(dat, transf);
      // pars.update_delta_eta1(dat, transf);
      // pars.update_delta_eta2(dat, transf);
      // pars.update_delta_eta_c(dat, transf);
      pars.update_beta(dat, transf);
      pars.update_delta_beta(dat, transf);
      pars.update_omega(dat, transf);
      pars.update_nu(dat, transf);
      pars.update_sigmasqeta(dat);
      // pars.update_a1234(dat);
    }
    progress_bar.increment();
    write_samples();
  }
  stop:
    NULL;
}

void SamplerWeak::write_samples() {
  pars.lambda_container.slice(pars.current_iter) = pars.lambda;
  pars.beta_container.slice(pars.current_iter) = pars.beta;
  pars.delta_beta_container.slice(pars.current_iter) = pars.delta_beta;
  // pars.delta_eta1_container.slice(pars.current_iter) = pars.delta_eta1;
  // pars.delta_eta2_container.slice(pars.current_iter) = pars.delta_eta2;
  pars.omega_container.col(pars.current_iter) = pars.omega;
  pars.xi_eta_container.slice(pars.current_iter) = pars.xi_eta;
  pars.zeta_container.col(pars.current_iter) = pars.zeta;
  pars.eta_container.slice(pars.current_iter) = pars.eta;
  pars.phi_container.slice(pars.current_iter) = pars.phi;
  pars.sigmasqetai_container.slice(pars.current_iter) = pars.sigmasqetai;
  pars.nu_container(pars.current_iter) = pars.nu;
  // pars.a1_container(pars.current_iter) = pars.a1;
  // pars.a2_container(pars.current_iter) = pars.a2;
  // pars.a3_container(pars.current_iter) = pars.a3;
  // pars.a4_container(pars.current_iter) = pars.a4;
  pars.nu_container(pars.current_iter) = pars.nu;
  pars.sigmasqeta_container.slice(pars.current_iter) = pars.sigmasqeta;
  pars.current_iter++;
}

Rcpp::List SamplerWeak::get_samples() {
  return Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("delta_beta", pars.delta_beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("omega", pars.omega_container),
                            Rcpp::Named("xi_eta", pars.xi_eta_container),
                            Rcpp::Named("zeta", pars.zeta_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("sigmasqetai", pars.sigmasqetai_container),
                            Rcpp::Named("sigmasqeta", pars.sigmasqeta_container),
                            Rcpp::Named("fit", transf.fit));
}

arma::vec SamplerWeak::get_grad(arma::vec& q) {
  arma::uword dim = q.n_elem;
  arma::mat delta_eta1 = arma::reshape(q.rows(0, dat.cdim * dat.nreg - 1), dat.nreg, dat.cdim);
  arma::mat delta_eta2 = arma::reshape(q.rows(dat.cdim * dat.nreg, dim - 1), dat.ldim, dat.cdim).t();
  arma::mat delta_eta1_cumprod = arma::mat(dat.nreg, dat.cdim);
  arma::mat delta_eta2_cumprod = arma::mat(dat.cdim, dat.ldim);
  for (arma::uword i = 0; i < dat.cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;
  arma::mat small_var;
  arma::vec etavec, betavec, etamean;
  arma::vec grad_vec = arma::vec(dat.cdim * (dat.nreg + dat.ldim));
  double grad;
  arma::vec delta_eta1_copy;
  arma::rowvec delta_eta2_copy;
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      delta_eta1_copy = delta_eta1.col(c);
      delta_eta1_copy(r) = 1;
      small_var = arma::cumprod(delta_eta1_copy) * arma::cumprod(delta_eta2.row(c));
      grad = 0;
      for (arma::uword rp = r; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = 0; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          grad = grad + .5 * dat.nsub * 1 / full_var(rp, lp) * small_var(rp, lp) -
            .5 * arma::as_scalar((etavec - etamean).t() *
            (etavec - etamean)) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(delta_eta1.col(c).row(r));
      if (r == 0) {
        // grad = grad + 1 / delta * (pars.a1 - 1) - 1;
      } else {
        // grad = grad + 1 / delta * (pars.a2 - 1) - 1;
      }
      grad_vec.row(r + dat.nreg * c) = -grad;
    }
  }
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      delta_eta2_copy = delta_eta2.row(c);
      delta_eta2_copy.col(l) = 1;
      small_var = arma::cumprod(delta_eta1.col(c)) * arma::cumprod(delta_eta2_copy);
      grad = 0;
      for (arma::uword rp = 0; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = l; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          grad = grad + .5 * dat.nsub * 1 / full_var(rp, lp) * small_var(rp, lp)- 
            .5 * arma::as_scalar((etavec - etamean).t() *
            (etavec - etamean)) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(delta_eta2.row(c).col(l));
      if (l == 0) {
        grad = grad + 1 / delta * (pars.a1 - 1) - 1;
      } else {
        grad = grad + 1 / delta * (pars.a2 - 1) - 1;
      }
      grad_vec(delta_eta1.n_cols * dat.nreg + c * dat.ldim + l) = -grad;
    }
  }
  // Rcpp::Rcout << grad_vec << "\n";
  return(grad_vec);
}

double SamplerWeak::get_density(arma::vec& q) {
  arma::uword dim = q.n_elem;
  arma::mat delta_eta1 = arma::reshape(q.rows(0, dat.cdim * dat.nreg - 1), dat.nreg, dat.cdim);
  arma::mat delta_eta2 = arma::reshape(q.rows(dat.cdim * dat.nreg, dim - 1), dat.ldim, dat.cdim).t();
  arma::mat delta_eta1_cumprod = arma::mat(dat.nreg, dat.cdim);
  arma::mat delta_eta2_cumprod = arma::mat(dat.cdim, dat.ldim);
  for (arma::uword i = 0; i < dat.cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;  
  arma::vec etavec, betavec, etamean;
  arma::vec density_vec = arma::vec(dat.cdim * (dat.nreg + dat.ldim));
  arma::vec sd = arma::vec(dat.nsub);
  double density = 0;
  for (arma::uword r = 0; r < dat.nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, dat.nreg, (dat.nsub - 1) * dat.nreg + r);
    for (arma::uword l = 0; l < dat.ldim; l++) {
      etavec = arma::vec(pars.eta.col(l)).rows(r_ind);
      betavec = arma::vec(pars.beta.col(l)).rows(r * dat.designdim, (r + 1) * dat.designdim - 1);
      etamean = dat.design * betavec;
      sd = arma::pow(arma::mat(pars.xi_eta.col(l)).rows(r_ind) * full_var(r, l), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
    }
  }
  return(-density);
}
arma::mat SamplerWeak::leapfrog(arma::uword num_steps, double step_size) {
  arma::uword dim = pars.delta_eta1.n_elem + pars.delta_eta2.n_elem;
  arma::vec p = ::sqrt(.1) * arma::randn(dim);
  arma::vec q = arma::vec(dim);
  arma::mat pq = arma::mat(dim, 2);
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      q.row(c * dat.nreg + r) = pars.delta_eta1.col(c).row(r);
    }
  }
  for (arma::uword c = 0; c < dat.cdim; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      q.row(dat.cdim * dat.nreg + c * dat.ldim + l) = 
        pars.delta_eta2.row(c).col(l);
    }
  }
  for (arma::uword i = 0; i < num_steps; i++) {
    p = p - step_size / 2. * get_grad(q);
    q = q + step_size * p;
    p = p - step_size / 2. * get_grad(q);
  }
  pq.col(0) = -p;
  pq.col(1) = q;
  return(pq);
}

void SamplerWeak::update_delta_eta(arma::uword num_steps, double step_size) {
  arma::mat new_leap = leapfrog(num_steps, step_size);
  arma::vec new_p = new_leap.col(0);
  arma::vec new_q = new_leap.col(1);
  // Rcpp::Rcout << new_q << "\n";
  // Rcpp::Rcout << q - new_q << "\n";
  // Rcpp::Rcout << get_density(q) << "\n" << get_density(new_q) << "\n";
  // Rcpp::Rcout << dat.cdim << "\n";
  // arma::vec d1 = arma::vec(new_q.rows(0, dat.nreg * dat.cdim - 1));
  // arma::vec d2 = arma::vec(new_q.rows(dat.nreg * dat.cdim, q.n_elem - 1));
  
  // arma::mat d1 = arma::reshape(arma::vec(new_q.rows(0, dat.nreg * dat.cdim - 1)), dat.nreg, dat.cdim);
  // arma::mat d2 = arma::reshape(arma::vec(new_q.rows(dat.nreg * dat.cdim, q.n_elem - 1)), dat.ldim, dat.cdim).t();
  // d1.resize(dat.nreg, dat.cdim);
  // d2.resize(dat.ldim, dat.cdim).t();
  // Rcpp::Rcout <<d1 << "\n" << pars.delta_eta1 << "\n";
  // Rcpp::Rcout << get_density(q) << "\n" << get_density(new_q);
  double A = R::runif(0, 1);
  // Rcpp::Rcout << A << "\n";
  double H_old = .1 * .5 * arma::as_scalar(p.t() * p) -get_density(q);
  double H_new = .1 * .5 * arma::as_scalar(new_p.t() * new_p) -get_density(new_q);
  // Rcpp::Rcout << H_old << "\n" << H_new << "\n";
  Rcpp::Rcout << -get_density(q) << "\n" << -get_density(new_q) << "\n";
  // Rcpp::Rcout << .5 * arma::as_scalar(p.t() * p) << "\n" << .5 * arma::as_scalar(new_p.t() * new_p) << "\n";
  Rcpp::Rcout << std::exp(H_new - H_old) << "\n";
  // Rcpp::Rcout << get_density(q) << "\n";
  // Rcpp::Rcout << get_grad(q) << "\n";
  Rcpp::Rcout << "H_old: " << H_old << "\n" << "H_new: " << H_new << "\n";
  // Rcpp::Rcout << -get_density(new_q) << "\n";
  // Rcpp::Rcout << new_q << "\n";
  if (A < std::exp(H_new - H_old)) {
    Rcpp::Rcout << "accepted\n";
    p = new_p;
    q = new_q;
    pars.delta_eta1 = arma::reshape(q.rows(0, dat.nreg * dat.cdim - 1), dat.nreg, dat.cdim);
    pars.delta_eta2 = arma::reshape(q.rows(dat.nreg * dat.cdim, q.n_elem - 1), dat.cdim, dat.ldim);
  }
}
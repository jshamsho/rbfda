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
  arma::vec p = arma::randn(dim);
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
  pars.lambda_container.slice(current_iter) = pars.lambda;
  pars.beta_container.slice(current_iter) = pars.beta;
  pars.delta_beta_container.slice(current_iter) = pars.delta_beta;
  pars.delta_eta_container.slice(current_iter) = pars.delta_eta;
  pars.omega_container.col(current_iter) = pars.omega;
  pars.xi_eta_container.slice(current_iter) = pars.xi_eta;
  pars.zeta_container.col(current_iter) = pars.zeta;
  pars.eta_container.slice(current_iter) = pars.eta;
  pars.phi_container(current_iter) = pars.phi;
  pars.sigmasqetai_container.slice(current_iter) = pars.sigmasqetai;
  pars.delta_eta_container.slice(current_iter) = pars.delta_eta;
  pars.nu_container(current_iter) = pars.nu;
  pars.a1_container(current_iter) = pars.a1;
  pars.a2_container(current_iter) = pars.a2;
  pars.a3_container(current_iter) = pars.a3;
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
      transf.complete_response(dat, pars);
      pars.update_lambda(dat, transf);
      pars.update_zeta(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_eta(dat, transf);
      // pars.update_xi_eta(dat, transf);
      // pars.update_delta_eta1(dat, transf);
      // pars.update_delta_eta2(dat, transf);
      pars.update_beta(dat, transf);
      pars.update_delta_beta(dat, transf);
      pars.update_omega(dat, transf);
      // pars.update_nu(dat, transf);
      // pars.update_a12(dat);
    }
    progress_bar.increment();
    write_samples();
  }
  stop:
    NULL;
}

void SamplerWeak::write_samples() {
  pars.lambda_container.slice(current_iter) = pars.lambda;
  pars.beta_container.slice(current_iter) = pars.beta;
  pars.delta_beta_container.slice(current_iter) = pars.delta_beta;
  pars.delta_eta1_container.slice(current_iter) = pars.delta_eta1;
  pars.delta_eta2_container.slice(current_iter) = pars.delta_eta2;
  pars.omega_container.col(current_iter) = pars.omega;
  pars.xi_eta_container.slice(current_iter) = pars.xi_eta;
  pars.zeta_container.col(current_iter) = pars.zeta;
  pars.eta_container.slice(current_iter) = pars.eta;
  pars.phi_container.slice(current_iter) = pars.phi;
  pars.sigmasqetai_container.slice(current_iter) = pars.sigmasqetai;
  pars.nu_container(current_iter) = pars.nu;
  pars.a1_container(current_iter) = pars.a1;
  pars.a2_container(current_iter) = pars.a2;
  current_iter++;
}

Rcpp::List SamplerWeak::get_samples() {
  return Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("delta_beta", pars.delta_beta_container),
                            Rcpp::Named("delta_eta1", pars.delta_eta1_container),
                            Rcpp::Named("delta_eta2", pars.delta_eta2_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("omega", pars.omega_container),
                            Rcpp::Named("xi_eta", pars.xi_eta_container),
                            Rcpp::Named("zeta", pars.zeta_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("sigmasqetai", pars.sigmasqetai_container),
                            Rcpp::Named("nu", pars.nu_container),
                            Rcpp::Named("a1", pars.a1_container),
                            Rcpp::Named("a2", pars.a2_container),
                            Rcpp::Named("fit", transf.fit));
}

arma::vec SamplerWeak::get_grad() {
  arma::mat delta_eta1_cumprod = arma::mat(dat.nreg, dat.c);
  arma::mat delta_eta2_cumprod = arma::mat(dat.c, dat.ldim);
  for (arma::uword i = 0; i < dat.c; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(pars.delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(pars.delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;
  arma::mat small_var;
  arma::vec etavec, betavec, etamean;
  arma::vec grad_vec = arma::vec(dat.c * (dat.nreg + dat.ldim));
  double delta_prod, grad;
  arma::vec delta_eta1_copy;
  arma::rowvec delta_eta2_copy;
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      delta_eta1_copy = pars.delta_eta1.col(c);
      delta_eta1_copy(r) = 1;
      small_var = arma::cumprod(delta_eta1_copy) * arma::cumprod(pars.delta_eta2.row(c));
      grad = 0;
      for (arma::uword rp = r; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = 0; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          grad = grad - .5 * arma::as_scalar((etavec - etamean).t() *
            arma::diagmat(arma::mat(pars.xi_eta.col(lp)).rows(r_ind)) * 
            (etavec - etamean)) * small_var(rp, lp) + 
            .5 * dat.nsub * 1 / full_var(rp, lp) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(pars.delta_eta1.col(c).row(r));
      if (r == 0) {
        grad = grad + 1 / delta * (pars.a1 - 1) - 1;
      } else {
        grad = grad + 1 / delta * (pars.a2 - 1) - 1;
      }
      grad_vec.row(r + dat.nreg * c) = -grad;
    }
  }
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      delta_eta2_copy = pars.delta_eta2.row(c);
      delta_eta2_copy.col(l) = 1;
      small_var = arma::cumprod(pars.delta_eta1.col(c)) * arma::cumprod(delta_eta2_copy);
      grad = 0;
      for (arma::uword rp = 0; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = l; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          grad = grad - .5 * arma::as_scalar((etavec - etamean).t() *
            arma::diagmat(arma::mat(pars.xi_eta.col(lp)).rows(r_ind)) * 
            (etavec - etamean)) * small_var(rp, lp) +
            .5 * dat.nsub * 1 / full_var(rp, lp) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(pars.delta_eta2.row(c).col(l));
      if (l == 0) {
        grad = grad + 1 / delta * (pars.a1 - 1) - 1;
      } else {
        grad = grad + 1 / delta * (pars.a2 - 1) - 1;
      }
      grad_vec(pars.delta_eta1.n_cols * dat.nreg + c * dat.ldim + l) = -grad;
    }
  }
  return(grad_vec);
}

arma::vec SamplerWeak::get_density() {
  arma::mat delta_eta1_cumprod = arma::mat(dat.nreg, dat.c);
  arma::mat delta_eta2_cumprod = arma::mat(dat.c, dat.ldim);
  for (arma::uword i = 0; i < dat.c; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(pars.delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(pars.delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;  arma::mat small_var;
  arma::vec etavec, betavec, etamean;
  arma::vec density_vec = arma::vec(dat.c * (dat.nreg + dat.ldim));
  arma::vec sd = arma::vec(dat.nsub);
  double delta_prod, density;
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      density = 0;
      for (arma::uword rp = r; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = 0; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          sd = arma::pow(arma::mat(pars.xi_eta.col(lp)).rows(r_ind) * full_var(rp, lp), -.5);
          density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
        }
      }
      double delta = arma::as_scalar(pars.delta_eta1.col(c).row(r));
      if (r == 0) {
        density = density + (pars.a1 - 1) * log(delta)  - delta;
      } else {
        density = density + (pars.a2 - 1) * log(delta) - delta;
      }
      density_vec.row(r + dat.nreg * c) = -density;
    }
  }
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      density = 0;
      for (arma::uword rp = 0; rp < dat.nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, dat.nreg, (dat.nsub - 1) * dat.nreg + rp);
        for (arma::uword lp = l; lp < dat.ldim; lp++) {
          etavec = arma::vec(pars.eta.col(lp)).rows(r_ind);
          betavec = arma::vec(pars.beta.col(lp)).rows(rp * dat.designdim, (rp + 1) * dat.designdim - 1);
          etamean = dat.design * betavec;
          sd = arma::pow(arma::mat(pars.xi_eta.col(lp)).rows(r_ind) * full_var(rp, lp), -.5);
          density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
        }
      }
      double delta = arma::as_scalar(pars.delta_eta2.row(c).col(l));
      if (l == 0) {
        density = density + log(delta) * (pars.a1 - 1) - 
          delta;
      } else {
        density = density + log(delta) * (pars.a2 - 1) - 
          delta;
      }
      density_vec.row(pars.delta_eta1.n_cols * dat.nreg + c * dat.ldim + l) = -density;
    }
  }
  return(density_vec);
}
arma::mat SamplerWeak::leapfrog(arma::uword num_steps, double step_size) {
  arma::uword dim = pars.delta_eta1.n_elem + pars.delta_eta2.n_elem;
  arma::vec p = arma::randn(dim);
  arma::vec q = arma::vec(dim);
  arma::mat pq = arma::mat(dim, 2);
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword r = 0; r < dat.nreg; r++) {
      q.row(c * dat.nreg + r) = pars.delta_eta1.col(c).row(r);
    }
  }
  for (arma::uword c = 0; c < dat.c; c++) {
    for (arma::uword l = 0; l < dat.ldim; l++) {
      q.row(dat.c * dat.nreg + c * dat.ldim + l) = 
        pars.delta_eta2.row(c).col(l);
    }
  }
  for (arma::uword i = 0; i < num_steps; i++) {
    p = p - step_size / 2. * get_grad();
    q = q + step_size * p;
    pars.delta_eta1 = arma::reshape(q.rows(0, dat.c * dat.nreg - 1), dat.nreg, dat.c);
    pars.delta_eta2 = arma::reshape(q.rows(dat.c * dat.nreg, dim), dat.c, dat.ldim);
    p = p - step_size / 2. * get_grad();
  }
  pq.col(0) = -p;
  pq.col(1) = q;
  return(pq);
}

void SamplerWeak::update_delta_eta(arma::uword num_steps, double step_size) {
  
}
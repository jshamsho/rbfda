#include "utils.h"
static double const log2pi = std::log(2.0 * M_PI);

arma::uvec arma_mod(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

arma::vec bayesreg(arma::vec const &b, arma::mat const &Q) {
  arma::vec z, w, v;
  z = arma::randn<arma::vec>(b.n_elem);
  arma::mat chol = arma::chol(Q, "lower");
  w = arma::solve(arma::trimatl(chol), b, arma::solve_opts::fast);
  v = arma::solve(arma::trimatu(chol.t()), w + z, arma::solve_opts::fast);
  return v;
}

arma::vec bayesreg_orth(arma::vec const &b, arma::mat const &Q, arma::mat const &lin_constr) {
  arma::vec z, w, v;
  arma::mat chol, Cbar, Ctilde;
  z = arma::randn<arma::vec>(b.n_elem);
  chol = arma::chol(Q, "lower");
  w = arma::solve(arma::trimatl(chol), b, arma::solve_opts::fast);
  v = arma::solve(arma::trimatu(chol.t()), w + z, arma::solve_opts::fast);
  Cbar = arma::solve(arma::trimatl(chol), lin_constr.t(), arma::solve_opts::fast);
  Ctilde = arma::solve(arma::trimatu(chol.t()), Cbar, arma::solve_opts::fast);
  v = v - Ctilde * arma::solve(lin_constr * Ctilde, lin_constr * v);
  return v;
}


/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

// // [[Rcpp::export]]
// arma::vec delta_eta_log_density(arma::mat& delta_eta1, arma::mat& delta_eta2,
//                                 arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
//                                 double a1, double a2,
//                                 arma::uword nsub, arma::uword nreg, arma::uword ldim,
//                                 arma::uword designdim, arma::mat& design) {
//   
// }

// [[Rcpp::export]]
double get_delta_eta_density(arma::mat& delta_eta1, arma::mat& delta_eta2,
                         arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                         arma::mat& design) {
  if (arma::any(arma::vectorise(delta_eta1) <= 0) || 
      arma::any(arma::vectorise(delta_eta2) <= 0)) {
    return(arma::datum::inf);
  }
  arma::uword d = design.n_cols;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword nsub = eta.n_rows / nreg;
  arma::uword cdim = delta_eta1.n_cols;
  arma::mat delta_eta1_cumprod = arma::mat(nreg, cdim);
  arma::mat delta_eta2_cumprod = arma::mat(cdim, ldim);
  for (arma::uword i = 0; i < cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;  
  arma::vec etavec, betavec, etamean;
  arma::vec density_vec = arma::vec(cdim * (nreg + ldim));
  arma::vec sd = arma::vec(nsub);
  double density = 0;
  for (arma::uword r = 0; r < nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, nreg, (nsub - 1) * nreg + r);
    for (arma::uword l = 0; l < ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(r * d, (r + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * full_var(r, l), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
    }
  }
  for (arma::uword c = 0; c < cdim; c++) {
    for (arma::uword r = 0; r < nreg; r++) {
      density = density + R::dgamma(delta_eta1(r, c), 2, 1, true);
    }
    for (arma::uword l = 0; l < ldim; l++) {
      density = density + R::dgamma(delta_eta2(c, l), 2, 1, true);
    }
  }
  // density = density + R::dgamma(arma::as_scalar(delta_eta1.col(c).row(r)),)
  return(-density);
}

// [[Rcpp::export]]
Rcpp::List get_delta_eta_proposal(arma::mat& delta_eta1, arma::mat& delta_eta2) {
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword cdim = delta_eta1.n_cols;
  arma::mat delta_eta1_proposal = arma::mat(nreg, cdim);
  arma::mat delta_eta2_proposal = arma::mat(cdim, ldim);
  double prop;
  for (arma::uword c = 0; c < cdim; c++) {
    for (arma::uword r = 0; r < nreg; r++) {
      prop = -1;
      while(prop <= 0) {
        prop = delta_eta1(r, c) + .000005 * R::rnorm(0, 1);
      }
      delta_eta1_proposal(r, c) = prop;
    }
    for (arma::uword l = 0; l < ldim; l++) {
      prop = -1;
      while(prop <= 0) {
        prop = delta_eta2(c, l) + .25 * R::rnorm(0, 1);
      }
      delta_eta2_proposal(c, l) = prop;
    }
  }
  return(Rcpp::List::create(Rcpp::Named("delta_eta1_proposal", delta_eta1_proposal),
                            Rcpp::Named("delta_eta2_proposal", delta_eta2_proposal)));
}

// [[Rcpp::export]]
double get_delta_eta1_grad(arma::mat& delta_eta1, arma::mat& delta_eta2,
                          arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                          arma::mat& design, arma::uword r, arma::uword c) {
  
  arma::vec etavec, betavec, etamean;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword nsub = eta.n_rows / nreg;
  arma::uword cdim = delta_eta1.n_cols;
  arma::uword d = design.n_cols;
  arma::mat delta_eta1_cumprod = arma::mat(nreg, cdim);
  arma::mat delta_eta2_cumprod = arma::mat(cdim, ldim);
  arma::vec delta_eta1_copy;
  delta_eta1_copy = delta_eta1.col(c);
  delta_eta1_copy(r) = 1;
  for (arma::uword i = 0; i < cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;
  arma::mat small_var = arma::cumprod(delta_eta1_copy) * arma::cumprod(delta_eta2.row(c));
  double grad = 0;
  for (arma::uword rp = r; rp < nreg; rp++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
    for (arma::uword lp = 0; lp < ldim; lp++) {
      etavec = arma::vec(eta.col(lp)).rows(r_ind);
      betavec = arma::vec(beta.col(lp)).rows(rp * d, (rp + 1) * d - 1);
      etamean = design * betavec;
      grad = grad + .5 * nsub * 1 / full_var(rp, lp) * small_var(rp, lp) -
                    .5 * arma::as_scalar((etavec - etamean).t() * 
      arma::diagmat(arma::vec(xi_eta.col(lp)).rows(r_ind)) *
      (etavec - etamean)) * small_var(rp, lp);
    }
  }
  grad = grad + (2 - 1) / delta_eta1(r, c) - 1;
  // double delta = arma::as_scalar(delta_eta1.col(c).row(r));
  return(-grad);
  
}

// [[Rcpp::export]]
double get_delta_eta2_grad(arma::mat& delta_eta1, arma::mat& delta_eta2,
                           arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                           arma::mat& design, arma::uword c, arma::uword l) {
  
  arma::vec etavec, betavec, etamean;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword nsub = eta.n_rows / nreg;
  arma::uword cdim = delta_eta1.n_cols;
  arma::uword d = design.n_cols;
  arma::mat delta_eta1_cumprod = arma::mat(nreg, cdim);
  arma::mat delta_eta2_cumprod = arma::mat(cdim, ldim);
  arma::vec delta_eta2_copy;
  delta_eta2_copy = delta_eta2.row(c);
  delta_eta2_copy(l) = 1;
  for (arma::uword i = 0; i < cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;
  arma::mat small_var = arma::cumprod(delta_eta1.col(c)) * arma::cumprod(delta_eta2_copy);
  double grad = 0;
  for (arma::uword rp = 0; rp < nreg; rp++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
    for (arma::uword lp = l; lp < ldim; lp++) {
      etavec = arma::vec(eta.col(lp)).rows(r_ind);
      betavec = arma::vec(beta.col(lp)).rows(rp * d, (rp + 1) * d - 1);
      etamean = design * betavec;
      grad = grad + .5 * nsub * 1 / full_var(rp, lp) * small_var(rp, lp) -
                    .5 * arma::as_scalar((etavec - etamean).t() * 
      arma::diagmat(arma::vec(xi_eta.col(lp)).rows(r_ind)) *
      (etavec - etamean)) * small_var(rp, lp);
    }
  }
  grad = grad + (2 - 1) / delta_eta2(c, l) - 1;
  // double delta = arma::as_scalar(delta_eta1.col(c).row(r));
  return(-grad);
  
}

void compute_sigmasqeta_weak(arma::mat& delta_eta1,
                        arma::mat& delta_eta2,
                        arma::mat& sigmasqeta) {
  arma::uword cdim = delta_eta1.n_cols;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::mat delta_eta1_cumprod = arma::mat(nreg, cdim);
  arma::mat delta_eta2_cumprod = arma::mat(cdim, ldim);
  for (arma::uword i = 0; i < cdim; i++) {
    delta_eta1_cumprod.col(i) = arma::cumprod(delta_eta1.col(i));
    delta_eta2_cumprod.row(i) = arma::cumprod(delta_eta2.row(i));
  }
  sigmasqeta = delta_eta1_cumprod * delta_eta2_cumprod;
}

void compute_sigmasqeta_partial(arma::mat& delta_eta, arma::mat& sigmasqeta) {
  arma::rowvec init; 
  arma::uword ldim = delta_eta.n_cols;
  init = arma::cumprod(delta_eta.row(0));
  for (arma::uword l = 0; l < ldim; l++) {
    sigmasqeta.col(l) = arma::cumprod(delta_eta.col(l));
    if (l > 0) {
      sigmasqeta.col(l) = sigmasqeta.col(l) *
        init(l - 1);
    }
  }
}

void compute_sigmasqetai(arma::mat& sigmasqeta,
                              arma::mat& xi_eta,
                              arma::mat& sigmasqetai) {
  arma::uword nreg = sigmasqeta.n_rows;
  arma::uword ldim = sigmasqeta.n_cols;
  arma::uword nsub = xi_eta.n_rows / nreg;
  arma::uword idx;
  for (arma::uword i = 0; i < nsub; i++) {
    for (arma::uword l = 0; l < ldim; l++) {
      for (arma::uword r = 0; r < nreg; r++) {
        idx = i * nreg + r;
        sigmasqetai(idx, l) = xi_eta(idx, l) * sigmasqeta(r, l);
      }
    }
  }
}

// 
// arma::mat leapfrog(arma::vec q, arma::vec p, arma::vec grad,
//                    arma::uword steps, double step_size) {
//   p = p - step_size / 2. * delta_eta_diff()
// }

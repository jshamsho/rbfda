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
  for (arma::uword rp = 0; rp < nreg; rp++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
    for (arma::uword lp = 0; lp < ldim; lp++) {
      etavec = arma::vec(eta.col(lp)).rows(r_ind);
      betavec = arma::vec(beta.col(lp)).rows(rp * d, (rp + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(lp)).rows(r_ind) * full_var(rp, lp), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
      if (rp == 0 && lp == 0) {
        Rcpp::Rcout << arma::stddev(etavec) << "\n";
      }
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
  return(density);
}

double get_delta_eta1_density(arma::mat& delta_eta1, arma::mat& delta_eta2,
                              arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                              arma::mat& design, arma::uword r, arma::uword c) {
  if (arma::any(arma::vectorise(delta_eta1) <= 0)) {
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
  for (arma::uword rp = r; rp < nreg; rp++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
    for (arma::uword l = 0; l < ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(rp * d, (rp + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * full_var(rp, l), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
    }
  }
  density = density + R::dgamma(delta_eta1(r, c), 2, 1, true);
  return(density);
}

double get_delta_eta2_density(arma::mat& delta_eta1, arma::mat& delta_eta2,
                              arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                              arma::mat& design, arma::uword c, arma::uword l) {
  if (arma::any(arma::vectorise(delta_eta2) <= 0)) {
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
    for (arma::uword lp = l; lp < ldim; lp++) {
      etavec = arma::vec(eta.col(lp)).rows(r_ind);
      betavec = arma::vec(beta.col(lp)).rows(r * d, (r + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(lp)).rows(r_ind) * full_var(r, lp), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
    }
  }
  density = density + R::dgamma(delta_eta2(c, l), 2, 1, true);
  return(density);
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

// [[Rcpp::export]]
double compute_delta_eta_density_c(arma::mat& delta_eta1, arma::mat& delta_eta2,
                                   arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                                   arma::mat& design, arma::uword c, double a1,
                                   double a2, double a3, double a4) {
  double density = 0;
  arma::vec etavec, betavec, etamean;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword nsub = eta.n_rows / nreg;
  arma::uword cdim = delta_eta1.n_cols;
  arma::uword d = design.n_cols;
  arma::mat delta_eta1_cumprod = arma::mat(nreg, cdim);
  arma::mat delta_eta2_cumprod = arma::mat(cdim, ldim);
  for (arma::uword c = 0; c < cdim; c++) {
    delta_eta1_cumprod.col(c) = arma::cumprod(delta_eta1.col(c));
    delta_eta2_cumprod.row(c) = arma::cumprod(delta_eta2.row(c));
  }
  arma::mat full_var = delta_eta1_cumprod * delta_eta2_cumprod;
  arma::vec sd;
  for (arma::uword r = 0; r < nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, nreg, (nsub - 1) * nreg + r);
    for (arma::uword l = 0; l < ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(r * d, (r + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * full_var(r, l), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
      // Rcpp::Rcout << "r = " << r << "\nl = " << l << " " <<arma::accu(arma::log_normpdf(etavec, etamean, sd)) << "\n";
    }
  }
  for (arma::uword r = 0; r < nreg; r++) {
    if (r == 0) {
      density = density + R::dgamma(delta_eta1(r, c), a1, 1, true);
    } else {
      density = density + R::dgamma(delta_eta1(r, c), a2, 1, true);
    }
  }
  for (arma::uword l = 0; l < ldim; l++) {
    if (l == 0) {
      density = density + R::dgamma(delta_eta2(c, l), a3, 1, true);
    } else {
      density = density + R::dgamma(delta_eta2(c, l), a4, 1, true);
    }
  }
  return(density);
}

// [[Rcpp::export]]
double compute_delta_eta_density_c2(arma::mat& delta_eta1, arma::mat& delta_eta2,
                                   arma::vec& delta_eta3, arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                                   arma::mat& design, arma::uword c) {
  double density = 0;
  arma::vec etavec, betavec, etamean;
  arma::uword nreg = delta_eta1.n_rows;
  arma::uword ldim = delta_eta2.n_cols;
  arma::uword nsub = eta.n_rows / nreg;
  arma::uword cdim = delta_eta1.n_cols;
  arma::uword d = design.n_cols;
  arma::mat sigmasqeta = delta_eta1 * arma::diagmat(delta_eta3) * delta_eta2;
  arma::vec sd;
  // Rcpp::Rcout << arma::diagmat(delta_eta3) << "\n";
  double cut;
  for (arma::uword r = 0; r < nreg; r++) {
    arma::uvec r_ind = arma::regspace<arma::uvec>(r, nreg, (nsub - 1) * nreg + r);
    for (arma::uword l = 0; l < ldim; l++) {
      etavec = arma::vec(eta.col(l)).rows(r_ind);
      betavec = arma::vec(beta.col(l)).rows(r * d, (r + 1) * d - 1);
      etamean = design * betavec;
      sd = arma::pow(arma::mat(xi_eta.col(l)).rows(r_ind) * sigmasqeta(r, l), -.5);
      density = density + arma::accu(arma::log_normpdf(etavec, etamean, sd));
    }
  }
  // Rcpp::Rcout << density << "\n";
  for (arma::uword r = 0; r < nreg; r++) {
    if (r == 0) {
      cut = 0;
    } else {
      cut = delta_eta1(r - 1, c);
    }
    density = density + R::dgamma(delta_eta1(r, c), 1, 1, true) - 
      R::pgamma(cut, 1, 1, false, true);
  }
  for (arma::uword l = 0; l < ldim; l++) {
    if (l == 0) {
      cut = 0;
    }
    else {
      cut = delta_eta2(c, l - 1);
    }
    density = density + R::dgamma(delta_eta2(c, l), 1, 1, true) - 
      R::pgamma(cut, 1, 1, false, true);
    // Rcpp::Rcout << "l = " << l << ": " << log(R::dgamma(delta_eta2(c, l), 1, 1, false) / 
      // (1 - R::pgamma(cut, 1, 1, true, false))) << "\n";
  }
  if (c > 0) {
    cut = delta_eta3(c - 1);
    density = density + R::dgamma(delta_eta3(c), 1, 1, true) -
      R::pgamma(cut, 1, 1, false, true);
  }
  return(density);
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


void compute_sigmasqeta_weak2(arma::mat& delta_eta1,
                             arma::mat& delta_eta2,
                             arma::mat& sigmasqeta) {
  
  sigmasqeta = delta_eta1 * delta_eta2;
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

// [[Rcpp::export]]
arma::vec identify_delta1(arma::vec& input) {
  arma::vec output = arma::vec(input.n_elem);
  output(0) = input(0);
  input = arma::cumprod(input);
  input = arma::sort(input);
  if (input.n_elem > 1) {
    for (arma::uword i = 1; i < input.n_elem; i++) {
      output.row(i) = input.row(i) / input.row(i - 1);
    }
  }
  return(output);
}

// [[Rcpp::export]]
arma::rowvec identify_delta2(arma::rowvec& input) {
  arma::rowvec output = arma::rowvec(input.n_elem);
  output.col(0) = input.col(0);
  input = arma::cumprod(input);
  input = arma::sort(input);
  if (input.n_elem > 1) {
    for (arma::uword i = 1; i < input.n_elem; i++) {
      output.col(i) = input.col(i) / input.col(i - 1);
    }
  }
  return(output);
}

// [[Rcpp::export]]
double gam_trunc_left(double a, double b,  double cut){ 
  double u, pcut, y; 
  
  pcut = R::pgamma(cut,a,b,1,0);
  if(pcut>0.99){
    return(cut);
  } 
  u = R::runif(0.0,1.0); 
  u = pcut + (1-pcut)*u; 
  y = R::qgamma(u, a, b, 1, 0);
  return y; 
} 

// [[Rcpp::export]]
double gam_trunc_right(double a, double b, double cut) {
  double u, pcut, y;
  
  pcut = R::pgamma(cut, a, b, 0, 0);
  if (pcut > .99) {
    return(cut);
  }
  u = R::runif(0, 1);
  u = pcut + (1 - pcut) * u;
  y = R::qgamma(u, a, b, 0, 0);
  return(y);
}

// [[Rcpp::export]]
double normal_trunc_left(double a, double b, double cut) { 
  double u, pcut, y; 
  
  pcut = R::pnorm5(cut, a, b, 1, 0);
  if(pcut>0.99){
    return(cut);
  } 
  u = R::runif(0.0,1.0); 
  u = pcut + (1-pcut)*u; 
  y = R::qnorm5(u, a, b, 1, 0);
  return y; 
} 
// 
// [[Rcpp::export]]
arma::mat reshape_nreg(arma::mat eta, arma::uword nsub, arma::uword nreg) {
  arma::uword ldim = eta.n_cols;
  arma::uword col_ind;
  arma::uvec seqr = arma::uvec(nsub);
  arma::mat etamat(nsub, nreg * ldim);
  arma::mat etacor;
  for (arma::uword l = 0; l < ldim; l++) {
    for (arma::uword r = 0; r < nreg; r++) {
      col_ind = l * nreg + r;
      seqr = arma::regspace<arma::uvec>(r, nreg, r + nreg * (nsub - 1));
      etamat.col(col_ind) = arma::vec(eta.col(l)).rows(seqr);
    }
  }
  return etamat;
}

Rcpp::List calculate_waic_partial(Rcpp::List result) {
  Rcpp::List control = result["control"];
  Rcpp::List data = result["data"];
  Rcpp::List samples = result["samples"];
  arma::mat Y = data["response"];
  arma::mat B = data["basis"];
  arma::uword burnin = control["burnin"];
  arma::uword iterations = control["iterations"];
  arma::mat omega = samples["omega"];
  arma::cube lambda = samples["lambda"];
  arma::cube eta = samples["eta"];
  arma::field<arma::cube> phi = samples["phi"];
  arma::uword nt = B.n_rows;
  arma::uword nreg = omega.n_rows;
  arma::uword nsub = eta.slice(0).n_rows / nreg;
  arma::uword ldim = phi(0).n_slices;
  arma::mat prediction(nt, nreg);
  arma::vec loglik = arma::vec(nt * nsub * nreg);
  arma::running_stat_vec<arma::vec> running_loglik;
  arma::uword start, end, start_eta, end_eta, start_loglik, end_loglik;
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    arma::mat sdmat = arma::repmat(arma::pow(omega.col(iter).t(), -.5) , nt, 1);
    for (arma::uword i = 0; i < nsub; i++) {
      prediction.zeros();
      start = i * nt;
      end = (i + 1) * nt - 1;
      start_eta = i * nreg;
      end_eta = (i + 1) * nreg - 1;
      start_loglik = i * nt * nreg;
      end_loglik = (i + 1) * nt * nreg - 1;
      for (arma::uword l = 0; l < ldim; l++) {
        prediction = prediction +
          B * lambda.slice(iter).col(l) *
          eta.slice(iter).col(l).rows(start_eta, end_eta).t() *
          phi(iter).slice(l).t();
      }
      loglik.rows(start_loglik, end_loglik) = 
        arma::vectorise(arma::log_normpdf(Y.rows(start, end), prediction, sdmat));
    }
    running_loglik(loglik);
  }
  
  double lppd = arma::accu(running_loglik.mean());
  double pwaic = arma::accu(running_loglik.var());
  double waic = -2 * (lppd - pwaic);
  arma::vec waic_vec = -2 * (running_loglik.mean() - running_loglik.var());
  double waic_se = std::sqrt(arma::var(waic_vec) * loglik.n_elem);
  return(Rcpp::List::create(Rcpp::Named("waic", waic),
                            Rcpp::Named("pwaic", pwaic),
                            Rcpp::Named("lppd", lppd),
                            Rcpp::Named("waic_se", waic_se),
                            Rcpp::Named("waic_vec", waic_vec)));
}

Rcpp::List calculate_waic_weak(Rcpp::List result) {
  Rcpp::List control = result["control"];
  Rcpp::List data = result["data"];
  Rcpp::List samples = result["samples"];
  arma::mat Y = data["response"];
  arma::mat B = data["basis"];
  arma::uword burnin = control["burnin"];
  arma::uword iterations = control["iterations"];
  arma::mat omega = samples["omega"];
  arma::cube lambda = samples["lambda"];
  arma::cube eta = samples["eta"];
  arma::cube phi = samples["phi"];
  arma::uword nt = B.n_rows;
  arma::uword nreg = omega.n_rows;
  arma::uword nsub = eta.slice(0).n_rows / nreg;
  arma::uword ldim = lambda.slice(0).n_cols;
  arma::mat prediction(nt, nreg);
  arma::vec loglik = arma::vec(nt * nsub * nreg);
  arma::running_stat_vec<arma::vec> running_loglik;
  arma::uword start, end, start_eta, end_eta, start_loglik, end_loglik;
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    arma::mat sdmat = arma::repmat(arma::pow(omega.col(iter).t(), -.5) , nt, 1);
    for (arma::uword i = 0; i < nsub; i++) {
      prediction.zeros();
      start = i * nt;
      end = (i + 1) * nt - 1;
      start_eta = i * nreg;
      end_eta = (i + 1) * nreg - 1;
      start_loglik = i * nt * nreg;
      end_loglik = (i + 1) * nt * nreg - 1;
      prediction = prediction +
        B * lambda.slice(iter) *
        eta.slice(iter).rows(start_eta, end_eta).t() *
        phi.slice(iter).t();
      
      loglik.rows(start_loglik, end_loglik) = 
        arma::vectorise(arma::log_normpdf(Y.rows(start, end), prediction, sdmat));
    }
    running_loglik(loglik);
  }
  
  double lppd = arma::accu(running_loglik.mean());
  double pwaic = arma::accu(running_loglik.var());
  double waic = -2 * (lppd - pwaic);
  arma::vec waic_vec = -2 * (running_loglik.mean() - running_loglik.var());
  double waic_se = std::sqrt(arma::var(waic_vec) * loglik.n_elem);
  return(Rcpp::List::create(Rcpp::Named("waic", waic),
                            Rcpp::Named("pwaic", pwaic),
                            Rcpp::Named("lppd", lppd),
                            Rcpp::Named("waic_se", waic_se),
                            Rcpp::Named("waic_vec", waic_vec)));
}

// [[Rcpp::export]]
Rcpp::List compare_waic(Rcpp::List waic1, Rcpp::List waic2) {
  double elpd1 = waic1["waic"];
  double elpd2 = waic2["waic"];
  double diff = elpd1 - elpd2;
  arma::vec elpd_vec1 = waic1["waic_vec"];
  arma::vec elpd_vec2 = waic2["waic_vec"];
  arma::vec diff_vec = elpd_vec1 - elpd_vec2;
  double diff_se = std::sqrt(arma::var(diff_vec) * elpd_vec1.n_elem);
  return(Rcpp::List::create(Rcpp::Named("diff", diff),
                            Rcpp::Named("diff_se", diff_se)));
}

// [[Rcpp::export]]
Rcpp::List calculate_waic(Rcpp::List result) {
  Rcpp::List control = result["control"];
  std::string covstruct = control["covstruct"];
  Rcpp::List waic;
  if (covstruct == "partial") {
    waic = calculate_waic_partial(result);
  } else if (covstruct == "weak") {
    waic = calculate_waic_weak(result);
  }
  return waic;
}

// 
// arma::mat leapfrog(arma::vec q, arma::vec p, arma::vec grad,
//                    arma::uword steps, double step_size) {
//   p = p - step_size / 2. * delta_eta_diff()
// }

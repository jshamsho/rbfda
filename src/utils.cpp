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
arma::vec delta_eta_diff(arma::mat& delta_eta1, arma::mat& delta_eta2,
                    arma::mat& eta, arma::mat& beta, arma::mat& xi_eta,
                    double a1, double a2,
                    arma::uword nsub, arma::uword nreg, arma::uword ldim,
                    arma::uword designdim, arma::mat& design) {
  arma::mat full_var = delta_eta1 * delta_eta2;
  arma::mat small_var;
  arma::vec etavec, betavec, etamean;
  arma::vec grad_vec = arma::vec(delta_eta1.n_cols * (nreg + ldim));
  double delta_prod, grad;
  arma::vec delta_eta1_copy;
  arma::rowvec delta_eta2_copy;
  for (arma::uword r = 0; r < nreg; r++) {
    for (arma::uword c = 0; c < delta_eta1.n_cols; c++) {
      delta_eta1_copy = delta_eta1.col(c);
      delta_eta1_copy(r) = 1;
      small_var = delta_eta1_copy * delta_eta2.row(c);
      grad = 0;
      for (arma::uword rp = r; rp < nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
        for (arma::uword lp = 0; lp < ldim; lp++) {
          etavec = arma::vec(eta.col(lp)).rows(r_ind);
          betavec = arma::vec(beta.col(lp)).rows(rp * designdim, (rp + 1) * designdim - 1);
          etamean = design * betavec;
          grad = grad + .5 * arma::as_scalar((etavec - etamean).t() *
            arma::diagmat(arma::mat(xi_eta.col(lp)).rows(r_ind)) * 
            (etavec - etamean)) * (1 / full_var(rp, lp)) * small_var(rp, lp) - 
            .5 * nsub * 1 / full_var(rp, lp) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(delta_eta1.col(c).row(r));
      if (r == 0) {
        grad = grad + log(delta) * (a1 - 1) - 
          delta;
      } else {
        grad = grad + log(delta) * (a2 - 1) - 
          delta;
      }
      grad_vec.row(r + nreg * c) = -grad;
    }
  }
  for (arma::uword l = 0; l < ldim; l++) {
    for (arma::uword c = 0; c < delta_eta1.n_cols; c++) {
      delta_eta2_copy = delta_eta2.row(c);
      delta_eta2_copy.col(l) = 1;
      small_var = delta_eta1.col(c) * delta_eta2_copy;
      grad = 0;
      for (arma::uword rp = 0; rp < nreg; rp++) {
        arma::uvec r_ind = arma::regspace<arma::uvec>(rp, nreg, (nsub - 1) * nreg + rp);
        for (arma::uword lp = l; lp < ldim; lp++) {
          etavec = arma::vec(eta.col(lp)).rows(r_ind);
          betavec = arma::vec(beta.col(lp)).rows(rp * designdim, (rp + 1) * designdim - 1);
          etamean = design * betavec;
          grad = grad + .5 * arma::as_scalar((etavec - etamean).t() *
            arma::diagmat(arma::mat(xi_eta.col(lp)).rows(r_ind)) * 
            (etavec - etamean)) * (1 / full_var(rp, lp)) * small_var(rp, lp) -
            .5 * nsub * 1 / full_var(rp, lp) * small_var(rp, lp);
        }
      }
      double delta = arma::as_scalar(delta_eta2.row(c).col(l));
      if (l == 0) {
        grad = grad + log(delta) * (a1 - 1) - 
          delta;
      } else {
        grad = grad + log(delta) * (a2 - 1) - 
          delta;
      }
      grad_vec(delta_eta1.n_cols * nreg + c * ldim + l) = -grad;
    }
  }
  return(grad_vec);
}
// 
// arma::mat leapfrog(arma::vec q, arma::vec p, arma::vec grad,
//                    arma::uword steps, double step_size) {
//   p = p - step_size / 2. * delta_eta_diff()
// }

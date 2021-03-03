#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double get_test_stat(arma::mat eta) {
  arma::mat etacov;
  etacov = arma::cov(eta);
  arma::uvec lower_indices = arma::trimatl_ind(arma::size(etacov), -1);
  return arma::accu(arma::square(etacov.elem(lower_indices)));
}

// [[Rcpp::export]]
arma::mat reshape_nreg(arma::mat eta, arma::uword nsub, arma::uword nreg) {
  arma::uword ldim = eta.n_cols;
  arma::uword col_ind;
  arma::uvec seqr = arma::uvec(nsub);
  arma::mat etamat(nsub, nreg * ldim);
  arma::mat etacov;
  for (arma::uword l = 0; l < ldim; l++) {
    for (arma::uword r = 0; r < nreg; r++) {
      col_ind = l * nreg + r;
      seqr = arma::regspace<arma::uvec>(r, nreg, r + nreg * (nsub - 1));
      etamat.col(col_ind) = arma::vec(eta.col(l)).rows(seqr);
    }
  }
  return etamat;
}
// [[Rcpp::export]]
arma::mat postcheck(Rcpp::List mcmc) {
  Rcpp::List samples = mcmc["samples"];
  Rcpp::List data = mcmc["data"];
  Rcpp::List control = mcmc["control"];
  arma::uword nsub = data["nsub"];
  arma::uword nreg = data["nreg"];
  arma::uword ldim = data["ldim"];
  arma::uword burnin = control["burnin"];
  arma::uword iterations = control["iterations"];
  arma::cube eta = samples["eta"];
  arma::cube sigmaetai = samples["sigmasqetai"];
  arma::cube xi_eta = samples["xi_eta"];
  arma::cube delta_eta = samples["delta_eta"];
  arma::mat etamat, sigmamat;
  arma::mat etagen = arma::mat(nsub, nreg * ldim);
  arma::mat test_stat = arma::mat(iterations - burnin, 2, arma::fill::zeros);
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    etamat = reshape_nreg(eta.slice(iter), nsub, nreg);
    sigmamat = reshape_nreg(sigmaetai.slice(iter), nsub, nreg);
    for (arma::uword i = 0; i < nsub; i++) {
      etagen.row(i) = arma::mvnrnd(arma::zeros(nreg * ldim), arma::diagmat(1. / sigmamat.row(i))).t();
    }
    test_stat(iter - burnin, 0) = get_test_stat(etamat);
    test_stat(iter - burnin, 1) = get_test_stat(etagen);
  }
  return test_stat;
}
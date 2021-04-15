#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// Define gamma1 here
double gamma1(double a, double b) {
  double g;
  g = ::tgamma((b + 1.) / 2.) * ::tgamma(a / 2.) /
    (::sqrt(arma::datum::pi) * ::tgamma((a + b) / 2.));
  return(g);
}
// [[Rcpp::export]]
double get_pval(arma::mat eta, arma::uword nreg, arma::uword ldim) {
  double p = eta.n_cols;
  double nsub = eta.n_rows;
  arma::mat etacor = arma::cor(eta);
  double dp = arma::accu(arma::abs(arma::trimatu(etacor, 1)));
  double mu = p * (p - 1) / 2 * gamma1(nsub - 1, 1);
  double tausq = p * (p - 1) / 2 * (1 / (nsub - 1) - ::pow(gamma1(nsub - 1, 1), 2));
  double stat = (dp - mu) / ::sqrt(tausq);
  double pval = R::pnorm(stat, 0.0, 1.0, false, false);
  return pval;
  // for (arma::uword l = 0; l < ldim - 1; l++) {
    // for (arma::uword r = 0; r < nreg; r++) {
      // sumsq = sumsq + arma::accu(arma::square(arma::vec(etacor.col(nreg * l + r)).rows((l + 1) * nreg + 1, nreg * ldim - 1)));
      // dp += arma::sum(arma::abs(etacor.col(nreg * l + r)))
        
    // }
  // }
  // arma::uvec lower_indices = arma::trimatl_ind(arma::size(etacor), -1);
  // return arma::accu(arma::square(etacor.elem(lower_indices)));
}

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
// [[Rcpp::export]]
Rcpp::List postcheck(Rcpp::List mcmc, arma::uword refdist_samples = 10000) {
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
  arma::vec test_stat = arma::vec(iterations - burnin);
  arma::vec null_test_stat = arma::vec(refdist_samples);
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    etamat = reshape_nreg(eta.slice(iter), nsub, nreg);
    // test_stat(iter - burnin) = get_test_stat(etamat, nreg, ldim);
  }
  for (arma::uword iter = 0; iter < refdist_samples; iter++) {
    for (arma::uword i = 0; i < nsub; i++) {
      etagen.row(i) = arma::randn<arma::rowvec>(nreg * ldim);
    }
    // null_test_stat(iter) = get_test_stat(etagen, nreg, ldim);
  }
  return Rcpp::List::create(Rcpp::Named("test_stat", test_stat),
                            Rcpp::Named("null_test_stat", null_test_stat));
}
#include "utils.h"

arma::uvec arma_mod(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}
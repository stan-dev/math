#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/prob/beta_rng.hpp>
#include <stan/math/prim/fun/read_corr_L.hpp>

namespace stan {
namespace math {

template <class RNG>
inline Eigen::MatrixXd lkj_corr_cholesky_rng(size_t K, double eta, RNG& rng) {
  static const char* function = "lkj_corr_cholesky_rng";
  check_positive(function, "Shape parameter", eta);

  Eigen::ArrayXd CPCs((K * (K - 1)) / 2);
  double alpha = eta + 0.5 * (K - 1);
  unsigned int count = 0;
  for (size_t i = 0; i < (K - 1); i++) {
    alpha -= 0.5;
    for (size_t j = i + 1; j < K; j++) {
      CPCs(count) = 2.0 * beta_rng(alpha, alpha, rng) - 1.0;
      count++;
    }
  }
  return read_corr_L(CPCs, K);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_RNG_HPP

#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/prob/beta_rng.hpp>
#include <stan/math/prim/meta/include_summand.hpp>
#include <stan/math/prim/fun/factor_cov_matrix.hpp>
#include <stan/math/prim/fun/factor_U.hpp>
#include <stan/math/prim/fun/read_corr_L.hpp>
#include <stan/math/prim/fun/read_corr_matrix.hpp>
#include <stan/math/prim/fun/read_cov_L.hpp>
#include <stan/math/prim/fun/read_cov_matrix.hpp>
#include <stan/math/prim/fun/make_nu.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/positive_constrain.hpp>
#include <stan/math/prim/fun/positive_free.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
#include <stan/math/prim/fun/lb_free.hpp>
#include <stan/math/prim/fun/ub_constrain.hpp>
#include <stan/math/prim/fun/ub_free.hpp>
#include <stan/math/prim/fun/lub_constrain.hpp>
#include <stan/math/prim/fun/lub_free.hpp>
#include <stan/math/prim/fun/prob_constrain.hpp>
#include <stan/math/prim/fun/prob_free.hpp>
#include <stan/math/prim/fun/corr_constrain.hpp>
#include <stan/math/prim/fun/corr_free.hpp>
#include <stan/math/prim/fun/simplex_constrain.hpp>
#include <stan/math/prim/fun/simplex_free.hpp>
#include <stan/math/prim/fun/ordered_constrain.hpp>
#include <stan/math/prim/fun/ordered_free.hpp>
#include <stan/math/prim/fun/positive_ordered_constrain.hpp>
#include <stan/math/prim/fun/positive_ordered_free.hpp>
#include <stan/math/prim/fun/cholesky_factor_constrain.hpp>
#include <stan/math/prim/fun/cholesky_factor_free.hpp>
#include <stan/math/prim/fun/cholesky_corr_constrain.hpp>
#include <stan/math/prim/fun/cholesky_corr_free.hpp>
#include <stan/math/prim/fun/corr_matrix_constrain.hpp>
#include <stan/math/prim/fun/corr_matrix_free.hpp>
#include <stan/math/prim/fun/cov_matrix_constrain.hpp>
#include <stan/math/prim/fun/cov_matrix_free.hpp>
#include <stan/math/prim/fun/cov_matrix_constrain_lkj.hpp>
#include <stan/math/prim/fun/cov_matrix_free_lkj.hpp>

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

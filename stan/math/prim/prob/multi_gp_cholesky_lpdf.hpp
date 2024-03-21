#ifndef STAN_MATH_PRIM_PROB_MULTI_GP_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_GP_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/mdivide_left_tri_low.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of a multivariate Gaussian Process for the given y, w, and
 * a Cholesky factor L of the kernel matrix Sigma.
 * Sigma = LL', a square, semi-positive definite matrix.
 * y is a dxN matrix, where each column is a different observation and each
 * row is a different output dimension.  The Gaussian Process is assumed to
 * have a scaled kernel matrix with a different scale for each output dimension.
 * This distribution is equivalent to:
 *    for (i in 1:d) row(y, i) ~ multi_normal(0, (1/w[i])*LL').
 *
 * @tparam T_y type of scalar
 * @tparam T_covar type of kernel
 * @tparam T_w type of weight
 *
 * @param y A dxN matrix
 * @param L The Cholesky decomposition of a kernel matrix
 * @param w A d-dimensional vector of positive inverse scale parameters for each
 * output.
 * @return The log of the multivariate GP density.
 * @throw std::domain_error if Sigma is not square, not symmetric,
 * or not semi-positive definite.
 */
template <bool propto, typename T_y, typename T_covar, typename T_w,
          require_all_eigen_matrix_dynamic_t<T_y, T_covar>* = nullptr,
          require_eigen_col_vector_t<T_w>* = nullptr>
return_type_t<T_y, T_covar, T_w> multi_gp_cholesky_lpdf(const T_y& y,
                                                        const T_covar& L,
                                                        const T_w& w) {
  using T_lp = return_type_t<T_y, T_covar, T_w>;
  static constexpr const char* function = "multi_gp_cholesky_lpdf";
  check_size_match(function, "Size of random variable (rows y)", y.rows(),
                   "Size of kernel scales (w)", w.size());
  check_size_match(function, "Size of random variable", y.cols(),
                   "rows of covariance parameter", L.rows());

  const auto& y_ref = to_ref(y);
  check_finite(function, "Random variable", y_ref);
  const auto& L_ref = to_ref(L);
  check_cholesky_factor(function, "Cholesky decomposition of kernel matrix",
                        L_ref);
  const auto& w_ref = to_ref(w);
  check_positive_finite(function, "Kernel scales", w_ref);

  if (y.rows() == 0) {
    return 0;
  }

  T_lp lp(0);
  if (include_summand<propto>::value) {
    lp += NEG_LOG_SQRT_TWO_PI * y.size();
  }

  if (include_summand<propto, T_covar>::value) {
    lp -= sum(log(L_ref.diagonal())) * y.rows();
  }

  if (include_summand<propto, T_w>::value) {
    lp += 0.5 * y.cols() * sum(log(w_ref));
  }

  if (include_summand<propto, T_y, T_w, T_covar>::value) {
    T_lp sum_lp_vec(0);
    for (int i = 0; i < y.rows(); i++) {
      sum_lp_vec
          += w_ref.coeff(i)
             * dot_self(mdivide_left_tri_low(L_ref, y_ref.row(i).transpose()));
    }
    lp -= 0.5 * sum_lp_vec;
  }

  return lp;
}

template <typename T_y, typename T_covar, typename T_w>
inline return_type_t<T_y, T_covar, T_w> multi_gp_cholesky_lpdf(const T_y& y,
                                                               const T_covar& L,
                                                               const T_w& w) {
  return multi_gp_cholesky_lpdf<false>(y, L, w);
}

}  // namespace math
}  // namespace stan
#endif

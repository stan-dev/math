#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CHECK_ARGS_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CHECK_ARGS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <cstddef>

namespace stan {
namespace math {
namespace internal {

/**
 * Checks the exponentially modified normal arguments in one branch-free pass
 * over the elements: y is not nan, mu is finite, and sigma and lambda are
 * positive and finite. The elementwise checks only run when that pass fails,
 * so they throw the same errors in the same order.
 *
 * @tparam T_y type of the random variable values
 * @tparam T_mu type of the location parameter values
 * @tparam T_sigma type of the scale parameter values
 * @tparam T_lambda type of the inverse scale parameter values
 * @param function name of the calling function, for error messages
 * @param y random variable values
 * @param mu location parameter values
 * @param sigma scale parameter values
 * @param lambda inverse scale parameter values
 * @return whether any element of y is infinite
 * @throw std::domain_error if any argument fails its check
 */
template <typename T_y, typename T_mu, typename T_sigma, typename T_lambda>
inline bool exp_mod_normal_check_args(const char* function, const T_y& y,
                                      const T_mu& mu, const T_sigma& sigma,
                                      const T_lambda& lambda) {
  bool good = true;
  bool y_inf = false;
  // With an empty container a scalar argument would be indexed past the
  // empty one, so the elementwise checks handle that case.
  if (!size_zero(y, mu, sigma, lambda)) {
    scalar_seq_view<T_y> y_vec(y);
    scalar_seq_view<T_mu> mu_vec(mu);
    scalar_seq_view<T_sigma> sigma_vec(sigma);
    scalar_seq_view<T_lambda> lambda_vec(lambda);
    const std::size_t N = max_size(y, mu, sigma, lambda);
    for (std::size_t n = 0; n < N; ++n) {
      const double y_n = value_of_rec(y_vec[n]);
      const double mu_n = value_of_rec(mu_vec[n]);
      const double sigma_n = value_of_rec(sigma_vec[n]);
      const double lambda_n = value_of_rec(lambda_vec[n]);
      good &= !std::isnan(y_n) & std::isfinite(mu_n) & (sigma_n > 0)
              & std::isfinite(sigma_n) & (lambda_n > 0)
              & std::isfinite(lambda_n);
      y_inf |= std::fabs(y_n) == INFTY;
    }
  } else {
    good = false;
  }
  if (unlikely(!good)) {
    check_not_nan(function, "Random variable", y);
    check_finite(function, "Location parameter", mu);
    check_positive_finite(function, "Scale parameter", sigma);
    check_positive_finite(function, "Inv_scale parameter", lambda);
  }
  return y_inf;
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif

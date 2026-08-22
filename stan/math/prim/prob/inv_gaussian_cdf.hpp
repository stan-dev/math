#ifndef STAN_MATH_PRIM_PROB_INV_GAUSSIAN_CDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAUSSIAN_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/inv_gaussian_lcdf.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the inverse Gaussian cumulative distribution function for the
 * given random variable, mean and shape. Given containers of matching sizes,
 * returns the product of probabilities.
 *
 * <p>The CDF is computed in log space and exponentiated once at the end,
 * since the factor \f$e^{2\lambda/\mu}\f$ in
 * \f$\Phi(z_1) + e^{2\lambda/\mu}\Phi(-z_2)\f$ overflows a double above
 * \f$2\lambda/\mu = 710\f$. See <code>inv_gaussian_lcdf</code>.
 *
 * <p>A boundary observation short-circuits the whole return:
 * \f$y = 0\f$ and \f$y = \infty\f$ each take the result to their
 * limit. The partials are zero at both.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of mean parameter
 * @tparam T_shape type of shape parameter
 * @param y scalar or container of scalars
 * @param mu mean parameter
 * @param lambda shape parameter
 * @return the product of probabilities
 * @throw std::domain_error if the mean or the shape is not positive and
 * finite, if the random variable is negative, or if any argument is NaN.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_y, typename T_loc, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_shape> inv_gaussian_cdf(
    const T_y& y, const T_loc& mu, const T_shape& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "inv_gaussian_cdf";
  check_consistent_sizes(function, "Random variable", y, "Mean parameter", mu,
                         "Shape parameter", lambda);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Mean parameter", mu_val);
  check_positive_finite(function, "Shape parameter", lambda_val);

  if (size_zero(y, mu, lambda)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, lambda_ref);

  if (sum(promote_scalar<int>(y_val == 0))) {
    return ops_partials.build(0.0);
  }
  if (sum(promote_scalar<int>(y_val == INFTY))) {
    return ops_partials.build(1.0);
  }

  constexpr bool any_ad = is_any_autodiff_v<T_y, T_loc, T_shape>;

  const auto& inv_mu = to_ref_if<any_ad>(inv(mu_val));
  const auto& inv_y = to_ref_if<any_ad>(inv(y_val));
  const auto& sqrt_lambda_over_y = to_ref_if<any_ad>(sqrt(lambda_val * inv_y));
  const auto& y_over_mu = to_ref_if<any_ad>(y_val * inv_mu);
  const auto& z1 = to_ref(sqrt_lambda_over_y * (y_over_mu - 1.0));
  const auto& z2 = to_ref(sqrt_lambda_over_y * (y_over_mu + 1.0));

  const auto& log_upper = to_ref(internal::log_scaled_upper_term(z1, z2));
  const auto& lcdf_elt = to_ref(log_sum_exp(internal::log_Phi(z1), log_upper));

  T_partials_return cdf = exp(sum(lcdf_elt));

  if constexpr (any_ad) {
    // d cdf / d theta_n = cdf * d log cdf_n / d theta_n
    //
    // Each partial is masked separately: at y == inf the shape partial is
    // 0 / 0 on its own. The inner select covers elements whose log CDF has
    // saturated to -inf at an interior y.
    const auto& w_dens = to_ref(cdf * exp(internal::log_phi(z1) - lcdf_elt));
    const auto& w_upper = to_ref(cdf * exp(log_upper - lcdf_elt));
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = w_dens * sqrt_lambda_over_y * inv_y;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -2.0 * lambda_val * w_upper * square(inv_mu);
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<2>(ops_partials)
          = 2.0 * w_upper * inv_mu - w_dens * inv_y / sqrt_lambda_over_y;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif

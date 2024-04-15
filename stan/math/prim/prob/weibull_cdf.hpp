#ifndef STAN_MATH_PRIM_PROB_WEIBULL_CDF_HPP
#define STAN_MATH_PRIM_PROB_WEIBULL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Weibull cumulative distribution function for the given
 * location and scale. Given containers of matching sizes, returns the
 * product of probabilities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if y is negative, alpha sigma is nonpositive
 */
template <typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> weibull_cdf(const T_y& y,
                                                 const T_shape& alpha,
                                                 const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "weibull_cdf";

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, alpha, sigma)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, sigma_ref);
  if (any(value_of_rec(y_val) == 0)) {
    return ops_partials.build(0.0);
  }

  constexpr bool any_derivs = !is_constant_all<T_y, T_shape, T_scale>::value;
  const auto& log_y = to_ref_if<any_derivs>(log(y_val));
  const auto& log_sigma = to_ref_if<any_derivs>(log(sigma_val));
  const auto& log_y_div_sigma = to_ref_if<any_derivs>(log_y - log_sigma);
  const auto& log_pow_n = to_ref_if<any_derivs>(alpha_val * log_y_div_sigma);
  const auto& pow_n = to_ref_if<any_derivs>(exp(log_pow_n));
  const auto& log_cdf_n = to_ref_if<any_derivs>(log1m_exp(-pow_n));

  T_partials_return log_cdf = sum(log_cdf_n);

  if (any_derivs) {
    const auto& log_rep_deriv = to_ref(log_pow_n + log_cdf - log_cdf_n - pow_n);
    if (!is_constant_all<T_y, T_scale>::value) {
      const auto& log_deriv_y_sigma = to_ref(log_rep_deriv + log(alpha_val));
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials) = exp(log_deriv_y_sigma - log_y);
      }
      if (!is_constant_all<T_scale>::value) {
        partials<2>(ops_partials) = -exp(log_deriv_y_sigma - log_sigma);
      }
    }
    if (!is_constant_all<T_shape>::value) {
      partials<1>(ops_partials) = exp(log_rep_deriv) * log_y_div_sigma;
    }
  }
  return ops_partials.build(exp(log_cdf));
}

}  // namespace math
}  // namespace stan
#endif

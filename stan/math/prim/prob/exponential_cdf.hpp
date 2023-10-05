#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the exponential cumulative distribution function for
 * the given y and beta.
 *
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 * @tparam T_y type of scalar
 * @tparam T_inv_scale type of inverse scale
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
 */
template <typename T_y, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_inv_scale>* = nullptr>
return_type_t<T_y, T_inv_scale> exponential_cdf(const T_y& y,
                                                const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_beta_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static const char* function = "exponential_cdf";
  T_y_ref y_ref = y;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(y, beta)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, beta_ref);

  constexpr bool any_derivatives = !is_constant_all<T_y, T_inv_scale>::value;
  const auto& exp_val = to_ref_if<any_derivatives>(exp(-beta_val * y_val));
  const auto& one_m_exp = to_ref_if<any_derivatives>(1 - exp_val);

  T_partials_return cdf(1.0);
  if (is_vector<T_y>::value || is_vector<T_inv_scale>::value) {
    cdf = forward_as<T_partials_array>(one_m_exp).prod();
  } else {
    cdf = forward_as<T_partials_return>(one_m_exp);
  }

  if (any_derivatives) {
    const auto& rep_deriv = to_ref_if<(
        !is_constant_all<T_y>::value && !is_constant_all<T_inv_scale>::value)>(
        exp_val / one_m_exp * cdf);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = beta_val * rep_deriv;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      partials<1>(ops_partials) = y_val * rep_deriv;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif

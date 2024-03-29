#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
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
 * The log of an exponential density for y with the specified
 * inverse scale parameter.
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y
 &\sim&
 \mbox{\sf{Expon}}(\beta) \\
 \log (p (y \, |\, \beta) )
 &=&
 \log \left( \beta \exp^{-\beta y} \right) \\
 &=&
 \log (\beta) - \beta y \\
 & &
 \mathrm{where} \; y > 0
 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_inv_scale type of inverse scale
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <bool propto, typename T_y, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_inv_scale>* = nullptr>
return_type_t<T_y, T_inv_scale> exponential_lpdf(const T_y& y,
                                                 const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_beta_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "exponential_lpdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(y, beta)) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, beta_ref);

  T_partials_return logp(0.0);
  if (include_summand<propto, T_inv_scale>::value) {
    logp = sum(log(beta_val)) * max_size(y, beta) / math::size(beta);
  }
  if (include_summand<propto, T_y, T_inv_scale>::value) {
    logp -= sum(beta_val * y_val);
  }

  if (!is_constant_all<T_y>::value) {
    using beta_val_scalar = scalar_type_t<decltype(beta_val)>;
    using beta_val_array = Eigen::Array<beta_val_scalar, Eigen::Dynamic, 1>;
    if (is_vector<T_y>::value && !is_vector<T_inv_scale>::value) {
      partials<0>(ops_partials) = T_partials_array::Constant(
          math::size(y), -forward_as<beta_val_scalar>(beta_val));
    } else if (is_vector<T_inv_scale>::value) {
      partials<0>(ops_partials) = -forward_as<beta_val_array>(beta_val);
    } else {
      forward_as<internal::broadcast_array<T_partials_return>>(
          partials<0>(ops_partials))
          = -forward_as<beta_val_scalar>(beta_val);
    }
  }
  if (!is_constant_all<T_inv_scale>::value) {
    partials<1>(ops_partials) = inv(beta_val) - y_val;
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_inv_scale>
inline return_type_t<T_y, T_inv_scale> exponential_lpdf(
    const T_y& y, const T_inv_scale& beta) {
  return exponential_lpdf<false>(y, beta);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_PROB_UNIFORM_LPDF_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of a uniform density for the given
 * y, lower, and upper bound.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{U}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta)) &=& \log \left( \frac{1}{\beta-\alpha}
 \right) \\
 &=& \log (1) - \log (\beta - \alpha) \\
 &=& -\log (\beta - \alpha) \\
 & & \mathrm{ where } \; y \in [\alpha, \beta], \log(0) \; \mathrm{otherwise}
 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_low type of lower bound
 * @tparam T_high type of upper bound
 * @param y A scalar variable.
 * @param alpha Lower bound.
 * @param beta Upper bound.
 * @throw std::invalid_argument if the lower bound is greater than
 *    or equal to the lower bound
 */
template <bool propto, typename T_y, typename T_low, typename T_high,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_low, T_high>* = nullptr>
return_type_t<T_y, T_low, T_high> uniform_lpdf(const T_y& y, const T_low& alpha,
                                               const T_high& beta) {
  using T_partials_return = partials_return_t<T_y, T_low, T_high>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_low>;
  using T_beta_ref = ref_type_if_not_constant_t<T_high>;
  static constexpr const char* function = "uniform_lpdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Lower bound parameter", alpha,
                         "Upper bound parameter", beta);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Lower bound parameter", alpha_val);
  check_finite(function, "Upper bound parameter", beta_val);
  check_greater(function, "Upper bound parameter", beta_val, alpha_val);

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_low, T_high>::value) {
    return 0.0;
  }
  if (sum(promote_scalar<int>(y_val < alpha_val))
      || sum(promote_scalar<int>(beta_val < y_val))) {
    return LOG_ZERO;
  }

  T_partials_return logp = 0;
  size_t N = max_size(y, alpha, beta);
  if (include_summand<propto, T_low, T_high>::value) {
    logp -= sum(log(beta_val - alpha_val)) * N / max_size(alpha, beta);
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  if (!is_constant_all<T_low, T_high>::value) {
    auto inv_beta_minus_alpha = to_ref_if<(!is_constant_all<T_high>::value
                                           && !is_constant_all<T_low>::value)>(
        inv(beta_val - alpha_val));
    if (!is_constant_all<T_high>::value) {
      if (is_vector<T_y>::value && !is_vector<T_low>::value
          && !is_vector<T_high>::value) {
        partials<2>(ops_partials) = -inv_beta_minus_alpha * math::size(y);
      } else {
        partials<2>(ops_partials) = -inv_beta_minus_alpha;
      }
    }
    if (!is_constant_all<T_low>::value) {
      if (is_vector<T_y>::value && !is_vector<T_low>::value
          && !is_vector<T_high>::value) {
        partials<1>(ops_partials) = inv_beta_minus_alpha * math::size(y);
      } else {
        partials<1>(ops_partials) = std::move(inv_beta_minus_alpha);
      }
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_low, typename T_high>
inline return_type_t<T_y, T_low, T_high> uniform_lpdf(const T_y& y,
                                                      const T_low& alpha,
                                                      const T_high& beta) {
  return uniform_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif

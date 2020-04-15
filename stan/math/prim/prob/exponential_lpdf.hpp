#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
template <bool propto, typename T_y, typename T_inv_scale>
return_type_t<T_y, T_inv_scale> exponential_lpdf(const T_y& y,
                                                 const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;
  using std::log;
  static const char* function = "exponential_lpdf";
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);

  if (size_zero(y, beta)) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_y, T_inv_scale> ops_partials(y, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t size_beta = stan::math::size(beta);
  size_t N = max_size(y, beta);

  VectorBuilder<include_summand<propto, T_inv_scale>::value, T_partials_return,
                T_inv_scale>
      log_beta(size_beta);
  if (include_summand<propto, T_inv_scale>::value) {
    for (size_t i = 0; i < size_beta; i++) {
      log_beta[i] = log(value_of(beta_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return y_dbl = value_of(y_vec[n]);
    if (include_summand<propto, T_inv_scale>::value) {
      logp += log_beta[n];
    }
    if (include_summand<propto, T_y, T_inv_scale>::value) {
      logp -= beta_dbl * y_dbl;
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= beta_dbl;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[n] += inv(beta_dbl) - y_dbl;
    }
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

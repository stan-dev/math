#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
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
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 * @tparam T_y Type of scalar.
 * @tparam T_inv_scale Type of inverse scale.
 */
template <bool propto, typename T_y, typename T_inv_scale>
inline auto exponential_lpdf(const T_y& y, const T_inv_scale& beta) {
  using T_partials = partials_return_t<T_y, T_inv_scale>;
  T_partials logp(0.0);

  using std::log;

  static const char* function = "exponential_lpdf";
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Random variable", y,
                         "Inverse scale parameter", beta);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_inv_scale> beta_vec(beta);
  const size_t N = max_size(y, beta);
  operands_and_partials<T_y, T_inv_scale> ops_partials(y, beta);
  if (size_zero(y, beta)) {
    return ops_partials.build(logp);
  }

  VectorBuilder<include_summand<propto, T_inv_scale>::value, T_partials,
                T_inv_scale>
      log_beta(length(beta));
  for (size_t i = 0; i < length(beta); i++) {
    if (include_summand<propto, T_inv_scale>::value) {
      log_beta[i] = log(value_of(beta_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials beta_dbl = value_of(beta_vec[n]);
    const T_partials y_dbl = value_of(y_vec[n]);
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
      ops_partials.edge2_.partials_[n] += 1 / beta_dbl - y_dbl;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_inv_scale>
inline auto exponential_lpdf(const T_y& y, const T_inv_scale& beta) {
  return exponential_lpdf<false>(y, beta);
}

}  // namespace math
}  // namespace stan
#endif

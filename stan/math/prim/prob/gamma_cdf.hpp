#ifndef STAN_MATH_PRIM_PROB_GAMMA_CDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The cumulative distribution function of the gamma distribution.
 * It is computed from the log of the regularized lower
 * incomplete gamma function, so it stays finite where the probability
 * underflows.
 *
 * @tparam T_y type of y
 * @tparam T_shape type of alpha
 * @tparam T_inv_scale type of beta
 * @param y random variable
 * @param alpha shape parameter
 * @param beta inverse scale parameter
 * @return probability, or product of probabilities for containers
 * @throw std::domain_error if a parameter is out of its support
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_shape, typename T_inv_scale>
inline return_type_t<T_y, T_shape, T_inv_scale> gamma_cdf(
    const T_y& y, const T_shape& alpha, const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_shape_ref = ref_type_t<T_shape>;
  using T_inv_scale_ref = ref_type_t<T_inv_scale>;
  static constexpr const char* function = "gamma_cdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_shape_ref alpha_ref = alpha;
  T_inv_scale_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);
  check_nonnegative(function, "Random variable", y_ref);

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  T_partials_return log_P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_shape_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_inv_scale_ref> beta_vec(beta_ref);
  const size_t N = max_size(y, alpha, beta);

  // An element with probability 0 (z = 0) makes the
  // result 0; the gradients are treated as zero
  for (size_t n = 0; n < N; n++) {
    const double y_rec = value_of_rec(y_vec.val(n));
    const double beta_rec = value_of_rec(beta_vec.val(n));
    if (beta_rec * y_rec == 0) {
      return ops_partials.build(0.0);
    }
  }

  constexpr bool need_d_a = is_autodiff_v<T_shape>;
  constexpr bool need_d_z = is_any_autodiff_v<T_y, T_inv_scale>;
  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_val = y_vec.val(n);
    const T_partials_return alpha_val = alpha_vec.val(n);
    const T_partials_return beta_val = beta_vec.val(n);
    const T_partials_return z = beta_val * y_val;
    const double z_rec = value_of_rec(z);
    // An element with probability 1 contributes nothing
    if (z_rec == INFTY) {
      continue;
    }
    const auto tail
        = internal::log_gamma_tail<false, need_d_a, need_d_z>(alpha_val, z);
    log_P += tail.log_prob;
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n] += beta_val * tail.d_z;
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<1>(ops_partials)[n] += tail.d_a;
    }
    if constexpr (is_autodiff_v<T_inv_scale>) {
      partials<2>(ops_partials)[n] += y_val * tail.d_z;
    }
  }
  const T_partials_return P = exp(log_P);
  if constexpr (is_autodiff_v<T_y>) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_shape>) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_inv_scale>) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      partials<2>(ops_partials)[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LPDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the inverse Gaussian density for the specified scalar(s) given
 * the specified mean(s) and shape(s). y, mu, or lambda can each be either a
 * scalar or a vector. Any vector inputs must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the log
 * probabilities for each observation/mean/shape triple.
 *
 * <p>The density is parameterized by the mean \f$\mu\f$ and the shape
 * \f$\lambda\f$ as \f$\sqrt{\lambda / (2\pi y^3)}
 * \exp(-\lambda (y - \mu)^2 / (2 \mu^2 y))\f$, with support
 * \f$y \in (0, \infty)\f$ and variance \f$\mu^3 / \lambda\f$.
 *
 * <p>Both \f$y = 0\f$ and \f$y = \infty\f$ are accepted and return
 * \f$-\infty\f$; only a negative \f$y\f$ throws.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of mean parameter
 * @tparam T_shape type of shape parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) mean parameter(s) for the inverse Gaussian
 * distribution.
 * @param lambda (Sequence of) shape parameter(s) for the inverse Gaussian
 * distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the mean or the shape is not positive and
 * finite, if the random variable is negative, or if any argument is NaN.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_loc, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_shape> inv_gaussian_lpdf(
    const T_y& y, const T_loc& mu, const T_shape& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "inv_gaussian_lpdf";
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
    return 0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_shape>::value) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, lambda_ref);

  // Both boundaries of the support have zero density and LOG_ZERO absorbs the
  // sum, so the whole container short-circuits exactly.
  // The gradients are technically ill-defined, but treated as zero.
  if (sum(promote_scalar<int>((y_val == 0) || (y_val == INFTY)))) {
    return ops_partials.build(LOG_ZERO);
  }

  constexpr bool any_ad = is_any_autodiff_v<T_y, T_loc, T_shape>;

  const auto& inv_mu = to_ref_if<any_ad>(inv(mu_val));
  const auto& inv_y = to_ref_if<any_ad>(inv(y_val));
  const auto& y_m_mu = to_ref_if<any_ad>(y_val - mu_val);
  // (y - mu)^2 / (2 mu^2 y)
  const auto& half_sq_scaled = to_ref(0.5 * square(y_m_mu * inv_mu) * inv_y);

  size_t N = max_size(y, mu, lambda);
  T_partials_return logp = -sum(lambda_val * half_sq_scaled);
  if constexpr (include_summand<propto>::value) {
    logp += N * NEG_LOG_SQRT_TWO_PI;
  }
  if constexpr (include_summand<propto, T_shape>::value) {
    logp += 0.5 * sum(log(lambda_val)) * N / math::size(lambda);
  }
  if constexpr (include_summand<propto, T_y>::value) {
    logp -= 1.5 * sum(log(y_val)) * N / math::size(y);
  }

  if constexpr (any_ad) {
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)
          = -1.5 * inv_y + 0.5 * lambda_val * (square(inv_y) - square(inv_mu));
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = lambda_val * y_m_mu * inv_mu * square(inv_mu);
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<2>(ops_partials) = 0.5 * inv(lambda_val) - half_sq_scaled;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_shape>
inline return_type_t<T_y, T_loc, T_shape> inv_gaussian_lpdf(
    const T_y& y, const T_loc& mu, const T_shape& lambda) {
  return inv_gaussian_lpdf<false>(y, mu, lambda);
}

}  // namespace math
}  // namespace stan
#endif

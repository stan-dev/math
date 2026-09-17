#ifndef STAN_MATH_PRIM_PROB_NORMAL_STANDARDIZE_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_STANDARDIZE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_scalar_ternary.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/** Standardize before evaluating a normal-family probability. Subtract first
 * to preserve precision near the mean; divide first only when subtraction
 * overflows for finite inputs. Decide before building autodiff nodes.
 */
template <typename T_y, typename T_mu, typename T_sigma>
inline auto normal_standardize(const T_y& y, const T_mu& mu,
                               const T_sigma& sigma) {
  if constexpr ((is_eigen_v<T_y> || is_eigen_v<T_mu> || is_eigen_v<T_sigma>)
                && !is_std_vector_v<T_y> && !is_std_vector_v<T_mu>
                && !is_std_vector_v<T_sigma>
                && std::is_arithmetic<scalar_type_t<T_y>>::value
                && std::is_arithmetic<scalar_type_t<T_mu>>::value
                && std::is_arithmetic<scalar_type_t<T_sigma>>::value) {
    // Preserve Eigen's vectorized arithmetic on the ordinary path. Only
    // overflowing differences need the scalar fallback.
    using Array = Eigen::Array<double, Eigen::Dynamic, 1>;
    const auto& diff = to_ref(as_array_or_scalar(promote_scalar<double>(y))
                              - as_array_or_scalar(mu));
    bool overflow;
    if constexpr (is_eigen_v<decltype(diff)>) {
      overflow = diff.size() && std::isinf(diff.abs().maxCoeff());
    } else {
      overflow = std::isinf(diff);
    }
    if (!overflow) {
      return Array(diff / as_array_or_scalar(sigma));
    }
    return Array(apply_scalar_ternary(
        [](const auto& yi, const auto& mi, const auto& si) {
          return normal_standardize(yi, mi, si);
        },
        y, mu, sigma));
  } else if constexpr (is_container<T_y>::value || is_container<T_mu>::value
                       || is_container<T_sigma>::value) {
    return apply_scalar_ternary(
        [](const auto& yi, const auto& mi, const auto& si) {
          return normal_standardize(yi, mi, si);
        },
        y, mu, sigma);
  } else if constexpr (std::is_integral<T_y>::value) {
    // Promote before subtraction or division without adding autodiff nodes.
    return normal_standardize(static_cast<double>(y), mu, sigma);
  } else {
    using R = return_type_t<T_y, T_mu, T_sigma>;
    const double yv = value_of_rec(y);
    const double mv = value_of_rec(mu);
    if (std::isinf(yv - mv) && std::isfinite(yv) && std::isfinite(mv)) {
      return R(y / sigma - mu / sigma);
    }
    return R((y - mu) / sigma);
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif

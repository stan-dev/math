#ifndef STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP
#define STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <utility>

namespace stan {
namespace math {
namespace internal {

// The erfcx pieces (Cody 1969 rationals and the small-argument polynomial)
// come from erfcx.hpp.

/** Tail factor in inv_x_sq = 1/x^2:
 * erfcx(x) = (1 + inv_x_sq * correction) / (sqrt(pi) x).
 */
template <typename T>
inline auto std_normal_tail_correction(const T& inv_x_sq) {
  // The arrays ascend in 1/x^2; reversed, they are Horner order.
  const auto pq = erfcx_tail_polynomials<erfcx_tail_p.size()>(
      erfcx_tail_p.crbegin(), erfcx_tail_q.crbegin(), inv_x_sq);
  return (pq.first / pq.second) / INV_SQRT_PI;
}

/** Scalar log Phi(z) and its slope phi(z) / Phi(z).
 * For z < 0.46875 sqrt(2) both come from erfcx with no exp; larger z needs
 * one exp for the complement. Infinite z gives -inf/0 values and inf/0
 * slopes.
 * The lower tail keeps log1p, inv_x_sq = 2 (1/neg_z)^2 and
 * neg_z / (1 + tail_term) so nested autodiff neither overflows nor loses the
 * 2 / neg_z^3 third derivative; the z > 40 return keeps derivatives finite
 * when x^2 overflows.
 */
template <bool calc_grad, typename T, require_stan_scalar_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(const T& z_in) {
  using T_return = return_type_t<T>;
  const auto z = z_in;
  if (z > 40) {
    if constexpr (calc_grad) {
      return std::pair<T_return, T_return>{0, 0};
    } else {
      return T_return(0);
    }
  } else if (z <= -4 * SQRT_TWO) {
    const auto neg_z = -z;
    const auto inv_x_sq = 2 * square(inv(neg_z));
    const auto tail_term = inv_x_sq * std_normal_tail_correction(inv_x_sq);
    const auto value
        = -(0.5 * z) * z - HALF_LOG_TWO_PI - log(neg_z) + log1p(tail_term);
    if constexpr (calc_grad) {
      return std::pair<T_return, T_return>{value, neg_z / (1.0 + tail_term)};
    } else {
      return T_return(value);
    }
  }
  // Not abs(z): its autodiff tangent is 0 at z == 0, which the slope needs.
  const auto x = (z < 0 ? -z : z) * INV_SQRT_TWO;
  if (x < 0.46875) {
    // Phi(z) = exp(-z^2 / 2) erfcx(-z / sqrt(2)) / 2 for either sign.
    const auto erfcx_neg = erfcx_small(-z * INV_SQRT_TWO);
    const auto value = LOG_HALF - (0.5 * z) * z + log(erfcx_neg);
    if constexpr (calc_grad) {
      return std::pair<T_return, T_return>{value,
                                           SQRT_TWO_OVER_SQRT_PI / erfcx_neg};
    } else {
      return T_return(value);
    }
  }
  const auto erfcx_x = x >= 4 ? erfcx_cody_tail(x) : erfcx_cody_middle(x);
  if (z < 0) {
    const auto value = -(0.5 * z) * z + LOG_HALF + log(erfcx_x);
    if constexpr (calc_grad) {
      return std::pair<T_return, T_return>{value,
                                           SQRT_TWO_OVER_SQRT_PI / erfcx_x};
    } else {
      return T_return(value);
    }
  }
  const auto density = exp(-(0.5 * z) * z);
  const auto tail = 0.5 * density * erfcx_x;
  // Not log1m(tail): its domain check costs ~3 ns per element here.
  const auto value = log1p(-tail);
  if constexpr (calc_grad) {
    return std::pair<T_return, T_return>{
        value, INV_SQRT_TWO_PI * density / (1.0 - tail)};
  } else {
    return T_return(value);
  }
}

/** Elementwise values and slopes; only the values when the gradient is off. */
template <bool calc_grad, typename T, require_eigen_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(T&& z) {
  using T_return = return_type_t<scalar_type_t<T>>;
  using Array = Eigen::Array<T_return, Eigen::Dynamic, 1>;
  if constexpr (calc_grad) {
    Array values(std::forward<T>(z));
    Array slopes(values.size());
    const auto value_size = values.size();
    for (Eigen::Index i = 0; i < value_size; ++i) {
      auto [value, grad]
          = std_normal_lcdf_value_grad<true>(T_return(values[i]));
      values[i] = value;
      slopes[i] = grad;
    }
    return std::make_pair(std::move(values), std::move(slopes));
  } else {
    Array values(std::forward<T>(z));
    for (Eigen::Index i = 0; i < values.size(); ++i) {
      values[i] = std_normal_lcdf_value_grad<false>(values[i]);
    }
    return values;
  }
}

/** Log of the standard normal cdf, or of its complement when `reflect` is
 * set: the input is negated so no negated autodiff operands are built.
 */
template <bool reflect, typename T_y>
inline return_type_t<T_y> std_normal_lcdf_impl(const char* function, T_y&& y) {
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  constexpr double sign = reflect ? -1.0 : 1.0;
  T_y_ref y_ref = std::forward<T_y>(y);
  // to_ref evaluates an expression argument once for the check and kernel.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  check_not_nan(function, "Random variable", y_val);

  if (size_zero(y_ref)) {
    return 0;
  }
  if constexpr (is_autodiff_v<T_y>) {
    const auto [values, slopes]
        = std_normal_lcdf_value_grad<true>(sign * y_val);
    auto ops_partials = make_partials_propagator(y_ref);
    partials<0>(ops_partials) = sign * slopes;
    return ops_partials.build(sum(values));
  } else {
    const auto values = std_normal_lcdf_value_grad<false>(sign * y_val);
    return sum(values);
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif

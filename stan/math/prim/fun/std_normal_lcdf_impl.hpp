#ifndef STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP
#define STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_square.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <array>
#include <utility>

namespace stan {
namespace math {
namespace internal {

// Rational approximations from Cody (1969), Math. Comp. 23:631-637 (SPECFUN
// CALERF). The erfcx coefficients are erfcx.hpp's, lowest power first.

inline constexpr std::array<double, 5> std_normal_erf_small_num
    = {1.85777706184603153e-1, 3.16112374387056560, 1.13864154151050156e2,
       3.77485237685302021e2, 3.20937758913846947e3};
// Highest power first, after a leading coefficient of 1.
inline constexpr std::array<double, 4> std_normal_erf_small_den
    = {2.36012909523441209e1, 2.44024637934444173e2, 1.28261652607737228e3,
       2.84423683343917062e3};

/** erf(x) for |x| < 0.46875, given x_sq = x^2. */
template <typename T, typename T_x_sq>
inline auto std_normal_erf_small(const T& x, const T_x_sq& x_sq) {
  constexpr auto& num = std_normal_erf_small_num;
  constexpr auto& den = std_normal_erf_small_den;
  const auto numerator
      = (((num[0] * x_sq + num[1]) * x_sq + num[2]) * x_sq + num[3]) * x_sq
        + num[4];
  const auto denominator
      = (((x_sq + den[0]) * x_sq + den[1]) * x_sq + den[2]) * x_sq + den[3];
  return x * numerator / denominator;
}

/** Tail factor in inv_x_sq = 1/x^2:
 * erfcx(x) = (1 + inv_x_sq * correction) / (sqrt(pi) x).
 */
template <typename T>
inline auto std_normal_tail_correction(const T& inv_x_sq) {
  constexpr auto& num = erfcx_tail_p;
  constexpr auto& den = erfcx_tail_q;
  const auto numerator
      = ((((num[5] * inv_x_sq + num[4]) * inv_x_sq + num[3]) * inv_x_sq
          + num[2])
             * inv_x_sq
         + num[1])
            * inv_x_sq
        + num[0];
  const auto denominator
      = ((((den[5] * inv_x_sq + den[4]) * inv_x_sq + den[3]) * inv_x_sq
          + den[2])
             * inv_x_sq
         + den[1])
            * inv_x_sq
        + den[0];
  return (numerator / denominator) / INV_SQRT_PI;
}

/** erfcx(x) = exp(x^2) erfc(x) for 0.46875 <= x <= 4. */
template <typename T>
inline auto std_normal_erfcx_middle(const T& x) {
  constexpr auto& num = erfcx_middle_p;
  constexpr auto& den = erfcx_middle_q;
  const auto numerator_high
      = ((((num[8] * x + num[7]) * x + num[6]) * x + num[5]) * x + num[4]) * x
        + num[3];
  const auto numerator
      = ((numerator_high * x + num[2]) * x + num[1]) * x + num[0];
  // den[8] is 1.
  const auto denominator_high
      = ((((x + den[7]) * x + den[6]) * x + den[5]) * x + den[4]) * x + den[3];
  const auto denominator
      = ((denominator_high * x + den[2]) * x + den[1]) * x + den[0];
  return numerator / denominator;
}

/** erfcx(x) = exp(x^2) erfc(x) for x > 4, given inv_x_sq = 1/x^2. */
template <typename T, typename T_inv_x_sq>
inline auto std_normal_erfcx_tail(const T& x, const T_inv_x_sq& inv_x_sq) {
  return (1 + inv_x_sq * std_normal_tail_correction(inv_x_sq)) * INV_SQRT_PI
         / x;
}

/** Scalar log Phi(z) and its slope phi(z) / Phi(z).
 * For z < 0 both come from erfcx with no exp; z > 0 needs one exp for the
 * complement. Infinite z gives -inf/0 values and inf/0 slopes.
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
    const auto erf_x = std_normal_erf_small(x, square(x));
    const auto value = LOG_HALF + log1p(z < 0 ? -erf_x : erf_x);
    if constexpr (calc_grad) {
      return std::pair<T_return, T_return>{
          value, SQRT_TWO_OVER_SQRT_PI * exp(-square(x))
                     / (z < 0 ? 1 - erf_x : 1 + erf_x)};
    } else {
      return T_return(value);
    }
  }
  const auto erfcx_x = x > 4 ? std_normal_erfcx_tail(x, inv_square(x))
                             : std_normal_erfcx_middle(x);
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
  decltype(auto) y_val = as_value_column_array_or_scalar(y_ref);
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

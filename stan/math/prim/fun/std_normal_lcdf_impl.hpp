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

/** erf(x) for |x| < 0.46875, given x2 = x^2. */
template <typename T, typename T2>
inline auto std_normal_erf_small(const T& x, const T2& x2) {
  constexpr auto& p = std_normal_erf_small_num;
  constexpr auto& q = std_normal_erf_small_den;
  const auto numerator
      = (((p[0] * x2 + p[1]) * x2 + p[2]) * x2 + p[3]) * x2 + p[4];
  const auto denominator = (((x2 + q[0]) * x2 + q[1]) * x2 + q[2]) * x2 + q[3];
  return x * numerator / denominator;
}

/** Tail factor in r = 1/x^2: erfcx(x) = (1 + r * correction) / (sqrt(pi) x). */
template <typename T>
inline auto std_normal_tail_correction(const T& r) {
  constexpr auto& p = erfcx_tail_p;
  constexpr auto& q = erfcx_tail_q;
  const auto numerator
      = ((((p[5] * r + p[4]) * r + p[3]) * r + p[2]) * r + p[1]) * r + p[0];
  const auto denominator
      = ((((q[5] * r + q[4]) * r + q[3]) * r + q[2]) * r + q[1]) * r + q[0];
  return (numerator / denominator) / INV_SQRT_PI;
}

/** erfcx(x) = exp(x^2) erfc(x) for 0.46875 <= x <= 4. */
template <typename T>
inline auto std_normal_erfcx_middle(const T& x) {
  constexpr auto& p = erfcx_middle_p;
  constexpr auto& q = erfcx_middle_q;
  const auto numerator_high
      = ((((p[8] * x + p[7]) * x + p[6]) * x + p[5]) * x + p[4]) * x + p[3];
  const auto numerator = ((numerator_high * x + p[2]) * x + p[1]) * x + p[0];
  // q[8] is 1.
  const auto denominator_high
      = ((((x + q[7]) * x + q[6]) * x + q[5]) * x + q[4]) * x + q[3];
  const auto denominator
      = ((denominator_high * x + q[2]) * x + q[1]) * x + q[0];
  return numerator / denominator;
}

/** erfcx(x) = exp(x^2) erfc(x) for x > 4, given r = 1/x^2. */
template <typename T, typename TR>
inline auto std_normal_erfcx_tail(const T& x, const TR& r) {
  return (1 + r * std_normal_tail_correction(r)) * INV_SQRT_PI / x;
}

/** Scalar log Phi(z) and its slope phi(z) / Phi(z).
 * For z < 0 both come from erfcx with no exp; z > 0 needs one exp for the
 * complement. Infinite z gives -inf/0 values and inf/0 slopes.
 * The lower tail keeps log1p, r = 2 (1/a)^2 and a / (1 + r c) so nested
 * autodiff neither overflows nor loses the 2 / a^3 third derivative; the
 * z > 40 return keeps derivatives finite when x^2 overflows.
 */
template <bool calc_grad, typename T, require_stan_scalar_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(const T& z_in) {
  using R = return_type_t<T>;
  const auto z = z_in;
  if (z > 40) {
    if constexpr (calc_grad) {
      return std::pair<R, R>{0, 0};
    } else {
      return R(0);
    }
  } else if (z <= -4 * SQRT_TWO) {
    const auto a = -z;
    const auto r = 2 * square(inv(a));
    const auto rc = r * std_normal_tail_correction(r);
    const auto value = -(0.5 * z) * z - HALF_LOG_TWO_PI - log(a) + log1p(rc);
    if constexpr (calc_grad) {
      return std::pair<R, R>{value, a / (1.0 + rc)};
    } else {
      return R(value);
    }
  }
  // Not abs(z): its autodiff tangent is 0 at z == 0, which the slope needs.
  const auto x = (z < 0 ? -z : z) * INV_SQRT_TWO;
  if (x < 0.46875) {
    const auto e = std_normal_erf_small(x, square(x));
    const auto value = LOG_HALF + log1p(z < 0 ? -e : e);
    if constexpr (calc_grad) {
      return std::pair<R, R>{value, SQRT_TWO_OVER_SQRT_PI * exp(-square(x))
                                        / (z < 0 ? 1 - e : 1 + e)};
    } else {
      return R(value);
    }
  }
  const auto errcx = x > 4 ? std_normal_erfcx_tail(x, inv_square(x))
                           : std_normal_erfcx_middle(x);
  if (z < 0) {
    const auto value = -(0.5 * z) * z + LOG_HALF + log(errcx);
    if constexpr (calc_grad) {
      return std::pair<R, R>{value, SQRT_TWO_OVER_SQRT_PI / errcx};
    } else {
      return R(value);
    }
  }
  const auto density = exp(-(0.5 * z) * z);
  const auto tail = 0.5 * density * errcx;
  // Not log1m(tail): its domain check costs ~3 ns per element here.
  const auto value = log1p(-tail);
  if constexpr (calc_grad) {
    return std::pair<R, R>{value, INV_SQRT_TWO_PI * density / (1.0 - tail)};
  } else {
    return R(value);
  }
}

/** Elementwise values and slopes; only the values when the gradient is off. */
template <bool calc_grad, typename T, require_eigen_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(T&& z) {
  using R = return_type_t<scalar_type_t<T>>;
  using Array = Eigen::Array<R, Eigen::Dynamic, 1>;
  if constexpr (calc_grad) {
    Array values(std::forward<T>(z));
    Array slopes(values.size());
    for (Eigen::Index i = 0; i < values.size(); ++i) {
      const auto result = std_normal_lcdf_value_grad<true>(R(values[i]));
      values[i] = result.first;
      slopes[i] = result.second;
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
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  check_not_nan(function, "Random variable", y_val);

  if (size_zero(y_ref)) {
    return 0;
  }

  auto z = sign * y_val;
  if constexpr (is_autodiff_v<T_y>) {
    const auto [values, slopes]
        = std_normal_lcdf_value_grad<true>(std::move(z));
    auto ops_partials = make_partials_propagator(y_ref);
    partials<0>(ops_partials) = sign * slopes;
    return ops_partials.build(sum(values));
  } else {
    const auto values = std_normal_lcdf_value_grad<false>(std::move(z));
    return sum(values);
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif

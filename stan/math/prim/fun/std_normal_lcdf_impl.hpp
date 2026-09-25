#ifndef STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP
#define STAN_MATH_PRIM_FUN_STD_NORMAL_LCDF_IMPL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_square.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <utility>

namespace stan {
namespace math {
namespace internal {

// Rational approximations from Cody (1969), Math. Comp. 23:631-637 (SPECFUN
// CALERF). Polynomial arithmetic only: no libm erf or erfc on the hot path.

/** erf(x) for |x| < 0.46875. */
template <typename T>
inline T std_normal_erf_small(const T& x) {
  static constexpr double a[]
      = {3.16112374387056560, 1.13864154151050156e2, 3.77485237685302021e2,
         3.20937758913846947e3, 1.85777706184603153e-1};
  static constexpr double b[] = {2.36012909523441209e1, 2.44024637934444173e2,
                                 1.28261652607737228e3, 2.84423683343917062e3};
  const T x2 = square(x);
  T numerator = a[4] * x2;
  T denominator = x2;
  for (int i = 0; i < 3; ++i) {
    numerator = (numerator + a[i]) * x2;
    denominator = (denominator + b[i]) * x2;
  }
  return x * (numerator + a[3]) / (denominator + b[3]);
}

/** Tail factor in r = 1/x^2: erfcx(x) = (1 + r * correction) / (sqrt(pi) x). */
template <typename T>
inline T std_normal_tail_correction(const T& r) {
  static constexpr double p[]
      = {0.000658749161529837803157, 0.0160837851487422766278,
         0.125781726111229246204,    0.360344899949804439429,
         0.305326634961232344035,    0.0163153871373020978498};
  static constexpr double q[]
      = {-0.00233520497626869185443, -0.0605183413124413191178,
         -0.527905102951428412248,   -1.87295284992346047209,
         -2.56852019228982242072,    -1.0};
  T numerator = p[5] * r + p[4];
  T denominator = q[5] * r + q[4];
  for (int i = 3; i >= 0; --i) {
    numerator = numerator * r + p[i];
    denominator = denominator * r + q[i];
  }
  return (numerator / denominator) / INV_SQRT_PI;
}

/** erfcx(x) = exp(x^2) erfc(x) for x >= 0.46875. */
template <typename T>
inline T std_normal_erfcx(const T& x) {
  if (x > 4) {
    const T r = inv_square(x);
    return (1 + r * std_normal_tail_correction(r)) * INV_SQRT_PI / x;
  }
  static constexpr double c[]
      = {5.64188496988670089e-1, 8.88314979438837594,   6.61191906371416295e1,
         2.98635138197400131e2,  8.81952221241769090e2, 1.71204761263407058e3,
         2.05107837782607147e3,  1.23033935479799725e3, 2.15311535474403846e-8};
  static constexpr double d[]
      = {1.57449261107098347e1, 1.17693950891312499e2, 5.37181101862009858e2,
         1.62138957456669019e3, 3.29079923573345963e3, 4.36261909014324716e3,
         3.43936767414372164e3, 1.23033935480374942e3};
  T numerator = c[8] * x;
  T denominator = x;
  for (int i = 0; i < 7; ++i) {
    numerator = (numerator + c[i]) * x;
    denominator = (denominator + d[i]) * x;
  }
  return (numerator + c[7]) / (denominator + d[7]);
}

/** Scalar log Phi(z) and its slope phi(z) / Phi(z).
 * For z < 0 both come from erfcx with no exp; z > 0 needs one exp for the
 * complement. Infinite z gives -inf/0 values and inf/0 slopes.
 * The lower tail keeps log1p, r = 2 (1/a)^2 and a / (1 + r c) so nested
 * autodiff neither overflows nor loses the 2 / a^3 third derivative; the
 * z > 40 return keeps derivatives finite when x^2 overflows.
 */
template <bool calc_grad, typename T, require_stan_scalar_t<T>* = nullptr>
inline std::pair<return_type_t<T>, return_type_t<T>> std_normal_lcdf_value_grad(
    const T& z_in) {
  using R = return_type_t<T>;
  const R z = z_in;
  if (z > 40) {
    return {0, 0};
  }
  if (z <= -4 * SQRT_TWO) {
    const R a = -z;
    const R r = 2 * square(inv(a));
    const R rc = r * std_normal_tail_correction(r);
    const R value = -(0.5 * z) * z - HALF_LOG_TWO_PI - log(a) + log1p(rc);
    if constexpr (calc_grad) {
      return {value, a / (1 + rc)};
    } else {
      return {value, 0};
    }
  }
  // Not abs(z): its autodiff tangent is 0 at z == 0, which the slope needs.
  const R x = (z < 0 ? R(-z) : z) * INV_SQRT_TWO;
  if (x < 0.46875) {
    const R e = std_normal_erf_small(x);
    const R value = LOG_HALF + log1p(z < 0 ? R(-e) : e);
    if constexpr (calc_grad) {
      return {value, SQRT_TWO_OVER_SQRT_PI * exp(-square(x))
                         / (z < 0 ? R(1 - e) : R(1 + e))};
    } else {
      return {value, 0};
    }
  }
  const R erfcx = std_normal_erfcx(x);
  if (z < 0) {
    const R value = -(0.5 * z) * z + LOG_HALF + log(erfcx);
    if constexpr (calc_grad) {
      return {value, SQRT_TWO_OVER_SQRT_PI / erfcx};
    } else {
      return {value, 0};
    }
  }
  const R density = exp(-(0.5 * z) * z);
  const R tail = 0.5 * density * erfcx;
  // Not log1m(tail): its domain check costs ~3 ns per element here.
  const R value = log1p(-tail);
  if constexpr (calc_grad) {
    return {value, INV_SQRT_TWO_PI * density / (1 - tail)};
  } else {
    return {value, 0};
  }
}

/** Elementwise values and slopes; empty slopes when the gradient is off. */
template <bool calc_grad, typename T, require_eigen_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(const T& z) {
  using R = return_type_t<scalar_type_t<T>>;
  using Array = Eigen::Array<R, Eigen::Dynamic, 1>;
  const auto& z_ref = to_ref(z);
  Array values(z_ref.size());
  Array slopes(calc_grad ? z_ref.size() : 0);
  for (Eigen::Index i = 0; i < z_ref.size(); ++i) {
    const auto result
        = std_normal_lcdf_value_grad<calc_grad>(R(z_ref.coeff(i)));
    values[i] = result.first;
    if constexpr (calc_grad) {
      slopes[i] = result.second;
    }
  }
  return std::make_pair(std::move(values), std::move(slopes));
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif

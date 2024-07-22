#ifndef STAN_MATH_PRIM_FUN_ASINH_HPP
#define STAN_MATH_PRIM_FUN_ASINH_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/arg.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/polar.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap `asinh()` so it can be vectorized.
 *
 * @tparam T argument scalar type
 * @param x argument
 * @return inverse hyperbolic sine of argument in radians.
 */
struct asinh_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::asinh;
    return asinh(x);
  }
};

/**
 * Returns the elementwise `asinh()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Inverse hyperbolic sine of each value in the container.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto asinh(const T& x) {
  return apply_scalar_unary<asinh_fun, T>::apply(x);
}

namespace internal {
/**
 * Return the hyperbolic arc sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc sine of the argument
 */
template <typename V>
inline std::complex<V> complex_asinh(const std::complex<V>& z) {
  std::complex<double> y_d = asinh(value_of_rec(z));
  auto y = log(z + sqrt(1 + z * z));
  return copysign(y, y_d);
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif

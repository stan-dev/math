#ifndef STAN_MATH_PRIM_FUN_ACOS_HPP
#define STAN_MATH_PRIM_FUN_ACOS_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/asin.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/prim/fun/polar.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap `acos()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arc cosine of variable in radians.
 */
struct acos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::acos;
    return acos(x);
  }
};

/**
 * Returns the elementwise `acos()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto acos(const Container& x) {
  return apply_scalar_unary<acos_fun, Container>::apply(x);
}

/**
 * Version of `acos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().acos(); });
}

namespace internal {
/**
 * Return the arc cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_acos(const std::complex<V>& z) {
  return 0.5 * pi() - asin(z);
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif

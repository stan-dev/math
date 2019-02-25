#ifndef STAN_MATH_PRIM_FUN_CBRT_HPP
#define STAN_MATH_PRIM_FUN_CBRT_HPP

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cube root of the specified value
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is negative.
 */
inline double cbrt(double x) { return std::cbrt(x); }

/**
 * Integer version of cbrt.
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double cbrt(int x) { return std::cbrt(x); }

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_CBRT_HPP
#define STAN_MATH_PRIM_FUN_CBRT_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/cbrt.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap cbrt() so it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Cube root of x.
 */
struct cbrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return cbrt(x);
  }
};

/**
 * Vectorized version of cbrt().
 * @param x Container of variables.
 * @tparam T Container type.
 * @return Cube root of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<cbrt_fun, T>::return_t cbrt(const T& x) {
  return apply_scalar_unary<cbrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_CBRT_HPP
#define STAN_MATH_PRIM_FUN_CBRT_HPP

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cube root of the specified value
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is negative.
 */
inline double cbrt(double x) { return std::cbrt(x); }

/**
 * Integer version of cbrt.
 *
 * @param[in] x Argument.
 * @return Cube root of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double cbrt(int x) { return std::cbrt(x); }

}  // namespace math
}  // namespace stan
#endif

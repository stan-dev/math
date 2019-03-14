#ifndef STAN_MATH_PRIM_FUN_LOG2_HPP
#define STAN_MATH_PRIM_FUN_LOG2_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>






namespace stan {
namespace math {

/**
 * Returns the base two logarithm of the argument (C99, C++11).
 *
 * The function is defined by:
 *
 * <code>log2(a) = log(a) / std::log(2.0)</code>.
 *
 * @param[in] u argument
 * @return base two logarithm of argument
 */
inline double log2(double u) {
  using std::log;
  return log(u) / LOG_2;
}

/**
 * Return the base two logarithm of the specified argument.  This
 * version is required to disambiguate <code>log2(int)</code>.
 *
 * @param[in] u argument
 * @return base two logarithm of argument
 */
inline double log2(int u) { return log2(static_cast<double>(u)); }

/**
 * Return natural logarithm of two.
 *
 * @return Natural logarithm of two.
 */
inline double log2() { return LOG_2; }














/**
 * Structure to wrap log2() so it can be vectorized.
 */
struct log2_fun {
  /**
   * Return the base two logarithm of the specified argument.
   *
   * @tparam T argument type
   * @param x argument
   * @return base two log of the argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    return log2(x);
  }
};

/**
 * Return the elementwise application of <code>log2()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T container type
 * @param x container
 * @return elementwise log2 of container elements
 */
template <typename T>
inline typename apply_scalar_unary<log2_fun, T>::return_t log2(const T& x) {
  return apply_scalar_unary<log2_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

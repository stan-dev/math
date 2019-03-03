#ifndef STAN_MATH_PRIM_FUN_LOG_HPP
#define STAN_MATH_PRIM_FUN_LOG_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/log.hpp>



namespace stan {
namespace math {

/**
 * Return the natural log of the specified argument.  This version
 * is required to disambiguate <code>log(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural log of argument.
 */
inline double log(int x) { return std::log(x); }

}  // namespace math
}  // namespace stan








namespace stan {
namespace math {

/**
 * Structure to wrap log() so that it can be vectorized.
 */
struct log_fun {
  /**
   * Return natural log of specified argument.
   *
   * @tparam T Scalar argument type.
   * @param[in] x Argument.
   * @return Natural log of x.
   */
  template <typename T>
  static inline T fun(const T& x) {
    using std::log;
    return log(x);
  }
};

/**
 * Return the elementwise natural log of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T Argument type.
 * @param[in] x Argument.
 * @return Elementwise application of natural log to the argument.
 */
template <typename T>
inline typename apply_scalar_unary<log_fun, T>::return_t log(const T& x) {
  return apply_scalar_unary<log_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_FUN_ATANH_HPP
#define STAN_MATH_PRIM_FUN_ATANH_HPP

#include <stanh/prim/fun/is_nan.hpp>
#include <stanh/prim/err/check_bounded.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic tangent of the specified value.
 * An argument of -1 returns negative infinity and an argument of 1
 * returns infinity.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is not in [-1, 1].
 */
inline double atanh(double x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_bounded("atanh", "x", x, -1.0, 1.0);
    return std::atanh(x);
  }
}

/**
 * Integer version of atanh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double atanh(int x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_bounded("atanh", "x", x, -1, 1);
    return std::atanh(x);
  }
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_ATANH_HPP
#define STAN_MATH_PRIM_FUN_ATANH_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/atanh.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap atanh() so it can be vectorized.
 */
struct atanh_fun {
  /**
   * Return the inverse hypberbolic tangent of the specified argument.
   *
   * @param x Argument.
   * @return Inverse hyperbolic tangent of the argument.
   * @tparam T Argument type.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return atanh(x);
  }
};

/**
 * Return the elementwise application of <code>atanh()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T Container type.
 * @param x Container.
 * @return Elementwise atanh of members of container.
 */
template <typename T>
inline typename apply_scalar_unary<atanh_fun, T>::return_t atanh(const T& x) {
  return apply_scalar_unary<atanh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_ATANH_HPP
#define STAN_MATH_PRIM_FUN_ATANH_HPP

#include <stanh/prim/fun/is_nan.hpp>
#include <stanh/prim/err/check_bounded.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic tangent of the specified value.
 * An argument of -1 returns negative infinity and an argument of 1
 * returns infinity.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is not in [-1, 1].
 */
inline double atanh(double x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_bounded("atanh", "x", x, -1.0, 1.0);
    return std::atanh(x);
  }
}

/**
 * Integer version of atanh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic tangent of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double atanh(int x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_bounded("atanh", "x", x, -1, 1);
    return std::atanh(x);
  }
}

}  // namespace math
}  // namespace stan
#endif

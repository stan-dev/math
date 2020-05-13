#ifndef STAN_MATH_PRIM_FUN_LAMBERT_W_HPP
#define STAN_MATH_PRIM_FUN_LAMBERT_W_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/lambert_w.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the W0 branch of the Lambert W function.
 *
 * @param[in] x argument
 * @return value of the W0 branch of the Lambert W function at argument
 * @throw std::domain_error if x is smaller than -e^(-1)
 */
inline double lambert_w0(double x) {
  return boost::math::lambert_w0(x, boost_policy_t());
}

/**
 * Return the value of the W0 branch of the Lambert W function.
 *
 * @param[in] x argument
 * @return value of the W0 branch of the Lambert W function at argument
 * @throw std::domain_error if x is smaller than -e^(-1)
 */
inline double lambert_w0(int x) {
  return boost::math::lambert_w0(x, boost_policy_t());
}

/**
 * Return the value of the W-1 branch of the Lambert W function.
 *
 * @param[in] x argument
 * @return value of the W-1 branch of the Lambert W function at argument
 * @throw std::domain_error if x is smaller than -e^(-1) or equal or bigger than 0
 */
inline double lambert_wm1(double x) {
  return boost::math::lambert_wm1(x, boost_policy_t());
}

/**
 * Return the value of the W-1 branch of the Lambert W function.
 *
 * @param[in] x argument
 * @return value of the W-1 branch of the Lambert W function at argument
 * @throw std::domain_error if x is smaller than -e^(-1) or equal or bigger than 0
 */
inline double lambert_wm1(int x) {
  return boost::math::lambert_wm1(x, boost_policy_t());
}

/**
 * Structure to wrap lambert_w0() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return value of the W0 branch of the Lambert W function at x.
 * @throw std::domain_error if x is smaller than -e^(-1)
 */
struct lambert_w0_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lambert_w0(x);
  }
};

/**
 * Vectorized version of lambert_w0().
 *
 * @tparam T type of container
 * @param x container
 * @return value of the W0 branch of the Lambert W function for each value in x
 * @throw std::domain_error if x is smaller than -e^(-1)
 */
template <typename T>
inline auto lambert_w0(const T& x) {
  return apply_scalar_unary<lambert_w0_fun, T>::apply(x);
}

/**
 * Structure to wrap lambert_wm1() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return value of the W-1 branch of the Lambert W function at x.
 * @throw std::domain_error if x is smaller than -e^(-1) or equal or bigger than 0
 */
struct lambert_wm1_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lambert_wm1(x);
  }
};

/**
 * Vectorized version of lambert_wm1().
 *
 * @tparam T type of container
 * @param x container
 * @return value of the W0 branch of the Lambert W function for each value in x
 * @throw std::domain_error if x is smaller than -e^(-1) or equal or bigger than 0
 */
template <typename T>
inline auto lambert_wm1(const T& x) {
  return apply_scalar_unary<lambert_wm1_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

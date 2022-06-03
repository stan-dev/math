#ifndef STAN_MATH_PRIM_FUN_LAMBERT_W_HPP
#define STAN_MATH_PRIM_FUN_LAMBERT_W_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/lambert_w.hpp>

namespace stan {
namespace math {

/**
 * Compute the Lambert W function on W0 branch for a value x.
 *
 * @tparam T type of value
 * @param x value
 * @return value of the W0 branch of the Lambert W function for x
 * @throw std::domain_error if x is less than or equal to `-e^(-1)`
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double lambert_w0(const T& x) {
  return boost::math::lambert_w0(x, boost_policy_t<>());
}

/**
 * Compute the Lambert W function on W-1 branch for a value x.
 *
 * @tparam T type of value
 * @param x value
 * @return value of the W-1 branch of the Lambert W function for x
 * @throw std::domain_error if x is less than or equal to `-e^(-1)` or greater
 * than or equal to 0
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double lambert_wm1(const T& x) {
  return boost::math::lambert_wm1(x, boost_policy_t<>());
}

namespace internal {

/**
 * Structure to wrap lambert_w0() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return value of the W0 branch of the Lambert W function at x.
 * @throw std::domain_error if x is less than or equal to `-e^(-1)`
 */
struct lambert_w0_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return lambert_w0(x);
  }
};

/**
 * Structure to wrap lambert_wm1() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return value of the W-1 branch of the Lambert W function at x.
 * @throw std::domain_error if x is less than or equal to `-e^(-1)` or greater
 * than or equal to 0
 */
struct lambert_wm1_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return lambert_wm1(x);
  }
};
}  // namespace internal

/**
 * Vectorized version of lambert_w0().
 *
 * @tparam T type of container
 * @param x container
 * @return value of the W0 branch of the Lambert W function for each value in x
 * @throw std::domain_error if x is less than or equal to `-e^(-1)`
 */
template <typename T, require_not_stan_scalar_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto lambert_w0(const T& x) {
  return apply_scalar_unary<internal::lambert_w0_fun, T>::apply(x);
}

/**
 * Vectorized version of lambert_wm1().
 *
 * @tparam T type of container
 * @param x container
 * @return value of the W0 branch of the Lambert W function for each value in x
 * @throw std::domain_error if x is less than or equal to `-e^(-1)` or greater
 * than or equal to 0
 */
template <typename T, require_not_stan_scalar_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto lambert_wm1(const T& x) {
  return apply_scalar_unary<internal::lambert_wm1_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_FUN_LOG1P_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of one plus the specified value.
 *
 * \f[
 * \mbox{log1p}(x) = \log(1 + x)
 * \f]
 *
 * This version is more stable for arguments near zero than
 * the direct definition.  If <code>x == -1</code>, <code>log1p(x)</code>
 * is defined to be negative infinity.
 *
 * @param[in] x Argument.
 * @return Natural log of one plus the argument.
 * @throw std::domain_error If argument is less than -1.
 */
inline double log1p(double x) {
  if (is_nan(x)) {
    return x;
  } else {
    check_greater_or_equal("log1p", "x", x, -1.0);
    return std::log1p(x);
  }
}

/**
 * Return the natural logarithm of one plus the specified
 * argument.  This version is required to disambiguate
 * <code>log1p(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural logarithm of one plus the argument.
 * @throw std::domain_error If argument is less than -1.
 */
inline double log1p(int x) {
  check_greater_or_equal("log1p", "x", x, -1);
  return std::log1p(x);
}

/**
 * Structure to wrap log1p() so it can be vectorized.
 */
struct log1p_fun {
  /**
   * Return the natural logarithm of one plus the specified value.
   *
   * @tparam T type of argument
   * @param x argument
   * @return natural log of one plus the argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    return log1p(x);
  }
};

/**
 * Return the elementwise application of <code>log1p()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise log1p of members of container.
 */
template <typename T,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline auto log1p(const T& x) {
  return apply_scalar_unary<log1p_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

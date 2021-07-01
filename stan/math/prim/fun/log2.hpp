#ifndef STAN_MATH_PRIM_FUN_LOG2_HPP
#define STAN_MATH_PRIM_FUN_LOG2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return natural logarithm of two.
 *
 * @return Natural logarithm of two.
 */
inline double log2() { return LOG_TWO; }

/**
 * Structure to wrap `log2()` so it can be vectorized.
 */
struct log2_fun {
  /**
   * Return the base two logarithm of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return base two log of the argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    using std::log2;
    return log2(x);
  }
};

/**
 * Return the elementwise application of `log2()` to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise log2 of container elements
 */
template <typename T, require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log2(const T& x) {
  return apply_scalar_unary<log2_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_LOG1M_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of one minus the specified value.
 *
 * The main use of this function is to cut down on intermediate
 * values during algorithmic differentiation.
 *
   \f[
   \mbox{log1m}(x) =
   \begin{cases}
     \ln(1-x) & \mbox{if } x \leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log1m}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{1-x} & \mbox{if } x \leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param[in] x Argument.
 * @return Natural log of one minus the argument.
 * @throw <code>std::domain_error</code> If the argument is greater than 1.
 * @throw <code>std::overflow_error</code> If the computation overflows.
 */
inline double log1m(double x) {
  if (!is_nan(x)) {
    check_less_or_equal("log1m", "x", x, 1);
  }
  return stan::math::log1p(-x);
}

/**
 * Structure to wrap log1m() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of (1 - x).
 */
struct log1m_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return log1m(x);
  }
};

/**
 * Vectorized version of log1m().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of 1 minus each value in x.
 */
template <
    typename T, require_not_var_matrix_t<T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log1m(const T& x) {
  return apply_scalar_unary<log1m_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

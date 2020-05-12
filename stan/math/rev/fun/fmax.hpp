#ifndef STAN_MATH_REV_FUN_FMAX_HPP
#define STAN_MATH_REV_FUN_FMAX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Returns the maximum of the two variable arguments (C99).
 *
 * No new variable implementations are created, with this function
 * defined as if by
 *
 * <code>fmax(a, b) = a</code> if a's value is greater than b's, and .
 *
 * <code>fmax(a, b) = b</code> if b's value is greater than or equal to a's.
 *
   \f[
   \mbox{fmax}(x, y) =
   \begin{cases}
     x & \mbox{if } x \geq y \\
     y & \mbox{if } x < y \\[6pt]
     x & \mbox{if } -\infty\leq x\leq \infty, y = \textrm{NaN}\\
     y & \mbox{if } -\infty\leq y\leq \infty, x = \textrm{NaN}\\
     \textrm{NaN} & \mbox{if } x, y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fmax}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } x \geq y \\
     0 & \mbox{if } x < y \\[6pt]
     1 & \mbox{if } -\infty\leq x\leq \infty, y = \textrm{NaN}\\
     0 & \mbox{if } -\infty\leq y\leq \infty, x = \textrm{NaN}\\
     \textrm{NaN} & \mbox{if } x, y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fmax}(x, y)}{\partial y} =
   \begin{cases}
     0 & \mbox{if } x \geq y \\
     1 & \mbox{if } x < y \\[6pt]
     0 & \mbox{if } -\infty\leq x\leq \infty, y = \textrm{NaN}\\
     1 & \mbox{if } -\infty\leq y\leq \infty, x = \textrm{NaN}\\
     \textrm{NaN} & \mbox{if } x, y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable.
 * @param b Second variable.
 * @return If the first variable's value is larger than the
 * second's, the first variable, otherwise the second variable.
 */
template <typename T>
inline var_value<T> fmax(const var_value<T>& a, const var_value<T>& b) {
  if (unlikely(is_nan(a))) {
    if (unlikely(is_nan(b))) {
      return var_value<T>(new precomp_vv_vari<T, vari_value<T>, vari_value<T>>(
          NOT_A_NUMBER, a.vi_, b.vi_, NOT_A_NUMBER, NOT_A_NUMBER));
    }
    return b;
  }
  if (unlikely(is_nan(b))) {
    return a;
  }
  return a > b ? a : b;
}

/**
 * Returns the maximum of the variable and scalar, promoting the
 * scalar to a variable if it is larger (C99).
 *
 * For <code>fmax(a, b)</code>, if a's value is greater than b,
 * then a is returned, otherwise a fresh variable implementation
 * wrapping the value b is returned.
 *
 * @param a First variable.
 * @param b Second value
 * @return If the first variable's value is larger than or equal
 * to the second value, the first variable, otherwise the second
 * value promoted to a fresh variable.
 */
template <typename T1, typename T2>
inline var_value<T1> fmax(const var_value<T1>& a, T2 b) {
  if (unlikely(is_nan(a))) {
    if (unlikely(is_nan(b))) {
      return var(new precomp_v_vari<T1, vari_value<T2>>(NOT_A_NUMBER, a.vi_,
                                                        NOT_A_NUMBER));
    }
    return var_value<T1>(b);
  }
  if (unlikely(is_nan(b))) {
    return a;
  }
  return a >= b ? a : var_value<T1>(b);
}

/**
 * Returns the maximum of a scalar and variable, promoting the scalar to
 * a variable if it is larger (C99).
 *
 * For <code>fmax(a, b)</code>, if a is greater than b's value,
 * then a fresh variable implementation wrapping a is returned, otherwise
 * b is returned.
 *
 * @param a First value.
 * @param b Second variable.
 * @return If the first value is larger than the second variable's value,
 * return the first value promoted to a variable, otherwise return the
 * second variable.
 */
template <typename T1, typename T2>
inline var_value<T2> fmax(T1 a, const var_value<T2>& b) {
  if (unlikely(is_nan(b))) {
    if (unlikely(is_nan(a))) {
      return var_value<T2>(new precomp_v_vari<T1, vari_value<T2>>(
          NOT_A_NUMBER, b.vi_, NOT_A_NUMBER));
    }
    return var_value<T2>(a);
  }
  if (unlikely(is_nan(a))) {
    return b;
  }
  return a > b ? var_value<T2>(a) : b;
}

}  // namespace math
}  // namespace stan
#endif

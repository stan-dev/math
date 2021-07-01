#ifndef STAN_MATH_PRIM_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_HPP

/**
 * The lgamma implementation in stan-math is based on either the
 * reentrant safe lgamma_r implementation from C or the
 * boost::math::lgamma implementation. The reentrant safe lgamma_r
 * implementation is preferred as it is faster when compared to the
 * boost version. However, the reentrant safe lgamma_r C function is
 * not available with MinGW compilers used on Windows. Therefore, the
 * boost version is used on Windows with MinGW compilers as fall
 * back. For details on the speed evaluations, please refer to
 * https://github.com/stan-dev/math/pull/1255 .
 */
#if !__MINGW32__
// _REENTRANT must be defined during compilation to ensure that cmath
// exports the reentrant safe lgamma_r version.
#if !_REENTRANT
#error \
    "stan-math requires _REENTRANT being defined during compilation" \
    "to make lgamma_r available."
#endif
#include <cmath>
#else
// MinGW compilers on Windows do not provide the reentrant lgamma_r
// such that we fall back to boost whenever we are on MinGW.
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <limits>
#endif
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the gamma function applied to
 * the specified argument.
 *
   \f[
   \mbox{lgamma}(x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \ln\Gamma(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lgamma}(x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \Psi(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
\f]
*
* @param x argument
* @return natural logarithm of the gamma function applied to
* argument
*/
inline double lgamma(double x) {
#if !__MINGW32__
  int sign = 1;
  return ::lgamma_r(x, &sign);
#else
  if (unlikely(x == 0.0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t<>());
#endif
}

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) {
#if !__MINGW32__
  int sign = 1;
  return ::lgamma_r(x, &sign);
#else
  if (unlikely(x == 0.0))
    return std::numeric_limits<double>::infinity();
  return boost::math::lgamma(x, boost_policy_t<>());
#endif
}

/**
 * Structure to wrap lgamma() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of the gamma function applied to x.
 * @throw std::domain_error if x is a negative integer or 0.
 */
struct lgamma_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lgamma(x);
  }
};

/**
 * Vectorized version of lgamma().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of the gamma function
 *         applied to each value in x.
 * @throw std::domain_error if any value is a negative integer or 0.
 */
template <typename T, require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto lgamma(const T& x) {
  return apply_scalar_unary<lgamma_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif

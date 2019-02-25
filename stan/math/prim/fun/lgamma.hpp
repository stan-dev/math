#ifndef STAN_MATH_PRIM_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_HPP

#include <cmath>

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
inline double lgamma(double x) { return std::lgamma(x); }

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) { return std::lgamma(x); }

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/lgamma.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap lgamma() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
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
 * @param x Container.
 * @tparam T Container type.
 * @return Natural log of the gamma function
 *         applied to each value in x.
 * @throw std::domain_error if any value is a negative integer or 0.
 */
template <typename T>
inline typename apply_scalar_unary<lgamma_fun, T>::return_t lgamma(const T& x) {
  return apply_scalar_unary<lgamma_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LGAMMA_HPP

#include <cmath>

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
inline double lgamma(double x) { return std::lgamma(x); }

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) { return std::lgamma(x); }

}  // namespace math
}  // namespace stan
#endif

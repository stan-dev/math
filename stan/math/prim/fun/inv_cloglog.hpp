#ifndef STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP

#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse complementary log-log function.
 *
 * The function is defined by
 *
 * <code>inv_cloglog(x) = 1 - exp(-exp(x))</code>.
 *
 * This function can be used to implement the inverse link
 * function for complementary-log-log regression.
 *
 *
   \f[
   \mbox{inv\_cloglog}(y) =
   \begin{cases}
     \mbox{cloglog}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv\_cloglog}(y)}{\partial y} =
   \begin{cases}
     \frac{\partial\, \mbox{cloglog}^{-1}(y)}{\partial y} & \mbox{if }
 -\infty\leq y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \mbox{cloglog}^{-1}(y) = 1 - \exp \left( - \exp(y) \right)
   \f]

   \f[
   \frac{\partial \, \mbox{cloglog}^{-1}(y)}{\partial y} = \exp(y-\exp(y))
   \f]
 *
 * @param x Argument.
 * @return Inverse complementary log-log of the argument.
 */
inline double inv_cloglog(double x) {
  using std::exp;
  return 1 - exp(-exp(x));
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP

#include <stanh/prim/vectorize/apply_scalar_unary.hpp>
#include <stanh/prim/fun/inv_cloglog.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_cloglog() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return 1 - exp(-exp(x)).
 */
struct inv_cloglog_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_cloglog(x);
  }
};

/**
 * Vectorized version of inv_cloglog().
 * @param x Container.
 * @tparam T Container type.
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_cloglog_fun, T>::return_t inv_cloglog(
    const T& x) {
  return apply_scalar_unary<inv_cloglog_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP

#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse complementary log-log function.
 *
 * The function is defined by
 *
 * <code>inv_cloglog(x) = 1 - exp(-exp(x))</code>.
 *
 * This function can be used to implement the inverse link
 * function for complementary-log-log regression.
 *
 *
   \f[
   \mbox{inv\_cloglog}(y) =
   \begin{cases}
     \mbox{cloglog}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv\_cloglog}(y)}{\partial y} =
   \begin{cases}
     \frac{\partial\, \mbox{cloglog}^{-1}(y)}{\partial y} & \mbox{if }
 -\infty\leq y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \mbox{cloglog}^{-1}(y) = 1 - \exp \left( - \exp(y) \right)
   \f]

   \f[
   \frac{\partial \, \mbox{cloglog}^{-1}(y)}{\partial y} = \exp(y-\exp(y))
   \f]
 *
 * @param x Argument.
 * @return Inverse complementary log-log of the argument.
 */
inline double inv_cloglog(double x) {
  using std::exp;
  return 1 - exp(-exp(x));
}

}  // namespace math
}  // namespace stan

#endif

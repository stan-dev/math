#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_CLOGLOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

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
 * @tparam T_x type of input kernel generator expression x
 * @param x Argument.
 * @return Inverse complementary log-log of the argument.
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline auto inv_cloglog(T_x&& x) {  // NOLINT
  return 1 - exp(-exp(std::forward<T_x>(x)));
}
}  // namespace math
}  // namespace stan

#endif
#endif

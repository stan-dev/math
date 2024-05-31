#ifndef STAN_MATH_PRIM_FUN_OWENS_T_HPP
#define STAN_MATH_PRIM_FUN_OWENS_T_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/owens_t.hpp>

namespace stan {
namespace math {

/**
 * Return the result of applying Owen's T function to the
 * specified arguments.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
   \f[
   \mbox{owens\_t}(h, a) =
   \begin{cases}
     \mbox{owens\_t}(h, a) & \mbox{if } -\infty\leq h, a \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } h = \textrm{NaN or } a = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{owens\_t}(h, a)}{\partial h} =
   \begin{cases}
     \frac{\partial\, \mbox{owens\_t}(h, a)}{\partial h} & \mbox{if }
 -\infty\leq h, a\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } h = \textrm{NaN
 or } a = \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{owens\_t}(h, a)}{\partial a} =
   \begin{cases}
     \frac{\partial\, \mbox{owens\_t}(h, a)}{\partial a} & \mbox{if }
 -\infty\leq h, a\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } h = \textrm{NaN
 or } a = \textrm{NaN} \end{cases} \f]

   \f[
   \mbox{owens\_t}(h, a) = \frac{1}{2\pi} \int_0^a
 \frac{\exp(-\frac{1}{2}h^2(1+x^2))}{1+x^2}dx \f]

   \f[
   \frac{\partial \, \mbox{owens\_t}(h, a)}{\partial h} =
 -\frac{1}{2\sqrt{2\pi}} \operatorname{erf}\left(\frac{ha}{\sqrt{2}}\right)
   \exp\left(-\frac{h^2}{2}\right)
   \f]

   \f[
   \frac{\partial \, \mbox{owens\_t}(h, a)}{\partial a} =
 \frac{\exp\left(-\frac{1}{2}h^2(1+a^2)\right)}{2\pi (1+a^2)} \f]
 *
 * @param h First argument
 * @param a Second argument
 * @return Owen's T function applied to the arguments.
 */
inline double owens_t(double h, double a) { return boost::math::owens_t(h, a); }

/**
 * Enables the vectorized application of the owens_t
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return owens_t function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_not_var_and_matrix_types<T1, T2>* = nullptr>
inline auto owens_t(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return owens_t(c, d); });
}

}  // namespace math
}  // namespace stan

#endif

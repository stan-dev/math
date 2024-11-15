#ifndef STAN_MATH_REV_FUN_CEIL_HPP
#define STAN_MATH_REV_FUN_CEIL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/ceil.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the ceiling of the specified variable (cmath).
 *
 * The derivative of the ceiling function is defined and
 * zero everywhere but at integers, and we set them to zero for
 * convenience,
 *
 * \f$\frac{d}{dx} {\lceil x \rceil} = 0\f$.
 *
 * The ceiling function rounds up.  For double values, this is the
 * smallest integral value that is not less than the specified
 * value.  Although this function is not differentiable because it
 * is discontinuous at integral values, its gradient is returned
 * as zero everywhere.
 *
   \f[
   \mbox{ceil}(x) =
   \begin{cases}
     \lceil x\rceil & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{ceil}(x)}{\partial x} =
   \begin{cases}
     0 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Input variable.
 * @return Ceiling of the variable.
 */
inline var ceil(const var& a) { return var(std::ceil(a.val())); }

template <typename T, require_matrix_t<T>* = nullptr>
inline auto ceil(const var_value<T>& a) {
  return var_value<T>(a.val().array().ceil());
}

}  // namespace math
}  // namespace stan
#endif

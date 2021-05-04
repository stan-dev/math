#ifndef STAN_MATH_FWD_FUNCTOR_APPLY_SCALAR_UNARY_HPP
#define STAN_MATH_FWD_FUNCTOR_APPLY_SCALAR_UNARY_HPP

#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/fwd/core/fvar.hpp>

namespace stan {
namespace math {

/**
 * Template specialization to fvar for vectorizing a unary scalar
 * function.  This is a base scalar specialization.  It applies
 * the function specified by the template parameter to the
 * argument.
 *
 * @tparam F Type of function to apply.
 * @tparam T Value and tangent type for for forward-mode
 * autodiff variable.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, fvar<T> > {
  /**
   * Function return type, which is same as the argument type for
   * the function, <code>fvar&lt;T&gt;</code>.
   */
  using return_t = fvar<T>;

  /**
   * Apply the function specified by F to the specified argument.
   *
   * @param x Argument variable.
   * @return Function applied to the variable.
   */
  static inline return_t apply(const fvar<T>& x) { return F::fun(x); }
};

}  // namespace math
}  // namespace stan
#endif

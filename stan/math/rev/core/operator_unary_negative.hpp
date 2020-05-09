#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_NEGATIVE_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_NEGATIVE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/v_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari>
class neg_vari : public op_vari<VariVal, Vari*> {
  using op_vari<VariVal, Vari*>::avi;
 public:
  explicit neg_vari(Vari* avi) : op_vari<VariVal, Vari*>(-(avi->val_), avi) {}
  void chain() {
    if (unlikely(is_nan(avi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ -= this->adj_;
    }
  }
};
}  // namespace internal

/**
 * Unary negation operator for variables (C++).
 *
 * \f$\frac{d}{dx} -x = -1\f$.
 *
   \f[
   \mbox{operator-}(x) =
   \begin{cases}
     -x & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator-}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Argument variable.
 * @return Negation of variable.
 */
template <typename T>
inline var_value<T> operator-(var_value<T> a) {
  return {new internal::neg_vari<T, vari_value<T>>(a.vi_)};
}

}  // namespace math
}  // namespace stan
#endif

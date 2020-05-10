#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_DECREMENT_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_DECREMENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari>
class decrement_vari : public op_vari<VariVal, Vari*> {
  using op_vari<VariVal, Vari*>::avi;

 public:
  explicit decrement_vari(Vari* avi)
      : op_vari<VariVal, Vari*>(avi->val_ - 1.0, avi) {}
  void chain() {
    if (unlikely(is_nan(avi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_;
    }
  }
};
}  // namespace internal

/**
 * Prefix decrement operator for variables (C++).
 *
 * Following C++, <code>(--a)</code> is defined to behave exactly as
 *
 * <code>a = a - 1.0)</code>
 *
 * does, but is faster and uses less memory.  In particular,
 * the result is an assignable lvalue.
 *
 * @param a Variable to decrement.
 * @return Reference the result of decrementing this input variable.
 */
template <typename T>
inline var_value<T>& operator--(var_value<T>& a) {
  a.vi_ = new internal::decrement_vari<T, vari_value<T>>(a.vi_);
  return a;
}

/**
 * Postfix decrement operator for variables (C++).
 *
 * Following C++, the expression <code>(a--)</code> is defined to
 * behave like the sequence of operations
 *
 * <code>var temp = a;  a = a - 1.0;  return temp;</code>
 *
 * @param a Variable to decrement.
 * @return Input variable.
 */
template <typename T>
inline var_value<T> operator--(var_value<T>& a, int /*dummy*/) {
  var_value<T> temp(a);
  a.vi_ = new internal::decrement_vari<T, vari_value<T>>(a.vi_);
  return temp;
}

}  // namespace math
}  // namespace stan
#endif

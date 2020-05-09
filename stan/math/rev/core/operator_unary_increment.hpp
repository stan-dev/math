#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_INCREMENT_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_INCREMENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/v_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari>
class increment_vari : public op_vari<VariVal, Vari*> {
  using op_vari<VariVal, Vari*>::avi;

 public:
  explicit increment_vari(Vari* avi)
      : op_vari<VariVal, Vari*>(avi->val_ + 1.0, avi) {}
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
 * Prefix increment operator for variables (C++).  Following C++,
 * (++a) is defined to behave exactly as (a = a + 1.0) does,
 * but is faster and uses less memory.  In particular, the
 * result is an assignable lvalue.
 *
 * @param a Variable to increment.
 * @return Reference the result of incrementing this input variable.
 */
template <typename T>
inline var_value<T>& operator++(var_value<T>& a) {
  a.vi_ = new internal::increment_vari<T, vari_value<T>>(a.vi_);
  return a;
}

/**
 * Postfix increment operator for variables (C++).
 *
 * Following C++, the expression <code>(a++)</code> is defined to behave like
 * the sequence of operations
 *
 * <code>var temp = a;  a = a + 1.0;  return temp;</code>
 *
 * @param a Variable to increment.
 * @return Input variable.
 */
template <typename T>
inline var_value<T> operator++(var_value<T>& a, int /*dummy*/) {
  var temp(a);
  a.vi_ = new internal::increment_vari<T, vari_value<T>>(a.vi_);
  return temp;
}

}  // namespace math
}  // namespace stan
#endif

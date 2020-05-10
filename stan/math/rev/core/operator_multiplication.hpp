#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class multiply_vari {};

template <typename VariVal, typename Vari1, typename Vari2>
class multiply_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    final : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  multiply_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ * bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += bvi()->val_ * this->adj_;
      bvi()->adj_ += avi()->val_ * this->adj_;
    }
  }
};

template <typename VariVal, typename Vari, typename Arith>
class multiply_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;

 public:
  multiply_vari(Vari* avi, const Arith& b)
      : op_vari<VariVal, Vari*, Arith>(avi->val_ * b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_ * bd();
    }
  }
};
}  // namespace internal

/**
 * Multiplication operator for two variables (C++).
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) = x\f$.
 *
   \f[
   \mbox{operator*}(x, y) =
   \begin{cases}
     xy & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator*}(x, y)}{\partial x} =
   \begin{cases}
     y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator*}(x, y)}{\partial y} =
   \begin{cases}
     x & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of multiplying operands.
 */
template <typename T>
inline var_value<T> operator*(const var_value<T>& a, const var_value<T>& b) {
  return {new internal::multiply_vari<T, vari_value<T>, vari_value<T>>(a.vi_,
                                                                       b.vi_)};
}

/**
 * Multiplication operator for a variable and a scalar (C++).
 *
 * The partial derivative for the variable is
 *
 * \f$\frac{\partial}{\partial x} (x * c) = c\f$, and
 *
 * @tparam Arith An arithmetic type
 * @param a Variable operand.
 * @param b Scalar operand.
 * @return Variable result of multiplying operands.
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline var_value<T> operator*(const var_value<T>& a, const Arith& b) {
  if (b == 1.0) {
    return a;
  }
  return {new internal::multiply_vari<T, vari_value<T>, Arith>(a.vi_, b)};
}

/**
 * Multiplication operator for a scalar and a variable (C++).
 *
 * The partial derivative for the variable is
 *
 * \f$\frac{\partial}{\partial y} (c * y) = c\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a Scalar operand.
 * @param b Variable operand.
 * @return Variable result of multiplying the operands.
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline var_value<T> operator*(const Arith& a, const var_value<T>& b) {
  if (a == 1.0) {
    return b;
  }
  return {new internal::multiply_vari<T, vari_value<T>, Arith>(
      b.vi_, a)};  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif

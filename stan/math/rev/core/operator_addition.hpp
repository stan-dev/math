#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/fill.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class add_vari {};


template <typename VariVal, typename Vari1, typename Vari2>
class add_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>> : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;
 public:
  add_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ + bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      fill(avi()->adj_, NOT_A_NUMBER);
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
      avi()->adj_ += this->adj_;
      bvi()->adj_ += this->adj_;
    }
  }
};

template <typename VariVal, typename Vari, typename Arith>
class add_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;
 public:
  add_vari(Vari* avi, Arith b) : op_vari<VariVal, Vari*, Arith>(avi->val_ + b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      fill(avi()->adj_, NOT_A_NUMBER);
    } else {
      avi()->adj_ += this->adj_;
    }
  }
};
}  // namespace internal

/**
 * Addition operator for variables (C++).
 *
 * The partial derivatives are defined by
 *
 * \f$\frac{\partial}{\partial x} (x+y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x+y) = 1\f$.
 *
 *
   \f[
   \mbox{operator+}(x, y) =
   \begin{cases}
     x+y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial y} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename T>
inline var_value<T> operator+(const var_value<T>& a, const var_value<T>& b) {
  return {new internal::add_vari<T, vari_value<T>, vari_value<T>>(a.vi_, b.vi_)};
}

/**
 * Addition operator for variable and scalar (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dx} (x + c) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First variable operand.
 * @param b Second scalar operand.
 * @return Result of adding variable and scalar.
 */
template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_value<T> operator+(const var_value<T>& a, Arith b) {
  if (b == 0.0) {
    return a;
  }
  return {new internal::add_vari<T, vari_value<T>, Arith>(a.vi_, b)};
}

/**
 * Addition operator for scalar and variable (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dy} (c + y) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First scalar operand.
 * @param b Second variable operand.
 * @return Result of adding variable and scalar.
 */
template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_value<T> operator+(Arith a, const var_value<T>& b) {
  if (a == 0.0) {
    return b;
  }
  return {new internal::add_vari<T, vari_value<T>, Arith>(b.vi_, a)};  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif

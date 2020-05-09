#ifndef STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/rev/core/dv_vari.hpp>
#include <stan/math/rev/core/op_vari.hpp>

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class subtract_vari {};

template <typename VariVal, typename Vari1, typename Vari2>
class subtract_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    final : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  subtract_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ - bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_;
      bvi()->adj_ -= this->adj_;
    }
  }
};

template <typename VariVal, typename Vari, typename Arith>
class subtract_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;

 public:
  subtract_vari(Vari* avi, const Arith& b)
      : op_vari<VariVal, Vari*, Arith>(avi->val_ - b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_;
    }
  }
};

template <typename VariVal, typename Arith, typename Vari>
class subtract_vari<
    VariVal, Arith, Vari,
    require_t<conjunction<std::is_arithmetic<Arith>, is_vari<Vari>>>>
    final : public op_vari<VariVal, Arith, Vari*> {
  using op_vari<VariVal, Arith, Vari*>::ad;
  using op_vari<VariVal, Arith, Vari*>::bvi;

 public:
  subtract_vari(const Arith& a, Vari* bvi)
      : op_vari<VariVal, Arith, Vari*>(a - bvi->val_, a, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(ad(), bvi()->val_))) {
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      bvi()->adj_ -= this->adj_;
    }
  }
};
}  // namespace internal

/**
 * Subtraction operator for variables (C++).
 *
 * The partial derivatives are defined by
 *
 * \f$\frac{\partial}{\partial x} (x-y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x-y) = -1\f$.
 *
   \f[
   \mbox{operator-}(x, y) =
   \begin{cases}
     x-y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator-}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator-}(x, y)}{\partial y} =
   \begin{cases}
     -1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam Var1 value type of a var
 * @tparam Var2 value type of a var
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting the second variable from
 * the first.
 */
template <typename T>
inline var_value<T> operator-(const var_value<T>& a, const var_value<T>& b) {
  return {new internal::subtract_vari<T, vari_value<T>, vari_value<T>>(a.vi_,
                                                                       b.vi_)};
}

/**
 * Subtraction operator for variable and scalar (C++).
 *
 * The derivative for the variable is
 *
 * \f$\frac{\partial}{\partial x} (x-c) = 1\f$, and
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param a First variable operand.
 * @param b Second scalar operand.
 * @return Result of subtracting the scalar from the variable.
 */
template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_value<T> operator-(const var_value<T>& a, const Arith& b) {
  if (b == 0.0) {
    return a;
  }
  return {new internal::subtract_vari<T, vari_value<T>, Arith>(a.vi_, b)};
}

/**
 * Subtraction operator for scalar and variable (C++).
 *
 * The derivative for the variable is
 *
 * \f$\frac{\partial}{\partial y} (c-y) = -1\f$, and
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param a First scalar operand.
 * @param b Second variable operand.
 * @return Result of subtracting a variable from a scalar.
 */
template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_value<T> operator-(const Arith& a, const var_value<T>& b) {
  return {new internal::subtract_vari<T, Arith, vari_value<T>>(a, b.vi_)};
}

}  // namespace math
}  // namespace stan
#endif

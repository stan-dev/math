#ifndef STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/rev/core/dv_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
class subtract_vv_vari final : public op_vv_vari {
 public:
  subtract_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(avi->val_ - bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
      bvi_->adj_ -= adj_;
    }
  }
};

class subtract_vd_vari final : public op_vd_vari {
 public:
  subtract_vd_vari(vari* avi, double b) : op_vd_vari(avi->val_ - b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
    }
  }
};

class subtract_dv_vari final : public op_dv_vari {
 public:
  subtract_dv_vari(double a, vari* bvi) : op_dv_vari(a - bvi->val_, a, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(ad_, bvi_->val_))) {
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      bvi_->adj_ -= adj_;
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
inline var operator-(const var& a, const var& b) {
  return {new internal::subtract_vv_vari(a.vi_, b.vi_)};
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator-(const var& a, Arith b) {
  if (b == 0.0) {
    return a;
  }
  return {new internal::subtract_vd_vari(a.vi_, b)};
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator-(Arith a, const var& b) {
  return {new internal::subtract_dv_vari(a, b.vi_)};
}

}  // namespace math
}  // namespace stan
#endif

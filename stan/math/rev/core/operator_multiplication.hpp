#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>

namespace stan {
namespace math {

namespace internal {
class multiply_vv_vari final : public op_vv_vari {
 public:
  multiply_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(avi->val_ * bvi->val_, avi, bvi) {}
  void chain() {
    avi_->adj_ += bvi_->val_ * adj_;
    bvi_->adj_ += avi_->val_ * adj_;
  }
};

class multiply_vd_vari final : public op_vd_vari {
 public:
  multiply_vd_vari(vari* avi, double b) : op_vd_vari(avi->val_ * b, avi, b) {}
  void chain() { avi_->adj_ += adj_ * bd_; }
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
inline var operator*(const var& a, const var& b) {
  return {new internal::multiply_vv_vari(a.vi_, b.vi_)};
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator*(const var& a, Arith b) {
  if (b == 1.0) {
    return a;
  }
  return {new internal::multiply_vd_vari(a.vi_, b)};
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator*(Arith a, const var& b) {
  if (a == 1.0) {
    return b;
  }
  return {new internal::multiply_vd_vari(b.vi_, a)};  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif

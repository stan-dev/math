#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

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
inline var operator+(const var& a, const var& b) {
  return make_callback_vari(a.vi_->val_ + b.vi_->val_,
                            [avi = a.vi_, bvi = b.vi_](const auto& vi) mutable {
                              if (unlikely(std::isnan(vi.val_))) {
                                avi->adj_ = NOT_A_NUMBER;
                                bvi->adj_ = NOT_A_NUMBER;
                              } else {
                                avi->adj_ += vi.adj_;
                                bvi->adj_ += vi.adj_;
                              }
                            });
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator+(const var& a, Arith b) {
  if (b == 0.0) {
    return a;
  }
  return make_callback_vari(a.vi_->val_ + b,
                            [avi = a.vi_, b](const auto& vi) mutable {
                              if (unlikely(std::isnan(vi.val_))) {
                                avi->adj_ = NOT_A_NUMBER;
                              } else {
                                avi->adj_ += vi.adj_;
                              }
                            });
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
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator+(Arith a, const var& b) {
  if (a == 0.0) {
    return b;
  }
  return b + a;  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif

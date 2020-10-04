#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fill.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/fun/add.hpp>


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
  var ret(a.val() + b.val());
  if (unlikely(is_any_nan(a.val(), b.val()))) {
    reverse_pass_callback([a, b]() mutable {
      a.adj() = NOT_A_NUMBER;
      b.adj() = NOT_A_NUMBER;
    });
  } else {
    reverse_pass_callback([ret, a, b]() mutable {
      a.adj() += ret.adj();
      b.adj() += ret.adj();
    });
  }
  return ret;
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
  } else {
    var ret(a.val() + b);
    if (unlikely(is_any_nan(a.val(), b))) {
      reverse_pass_callback([a]() mutable { a.adj() = NOT_A_NUMBER; });
    } else {
      reverse_pass_callback([ret, a]() mutable { a.adj() += ret.adj(); });
    }
    return ret;
  }
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
  return b + a;
}

/**
 * Addition operator for matrix variables (C++).
 *
 * @tparam VarMat1 A matrix of vars or a var with an underlying matrix type.
 * @tparam VarMat2 A matrix of vars or a var with an underlying matrix type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename VarMat1, typename VarMat2,
          require_any_var_matrix_t<VarMat1, VarMat2>* = nullptr>
inline auto operator+(const VarMat1& a, const VarMat2& b) {
  return add(a, b);
}

}  // namespace math
}  // namespace stan
#endif

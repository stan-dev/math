#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_NEGATIVE_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_NEGATIVE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

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
inline var operator-(const var& a) {
  return make_callback_var(
      -a.val(), [a](const auto& vi) mutable { a.adj() -= vi.adj(); });
}

/**
 * Compute additive inverse of input
 *
 * @tparam T type of input
 * @param a input
 * @return negative of input
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto operator-(const T& a) {
  return make_callback_var(-a.val(), [a](const auto& vi) mutable {
    for (Eigen::Index j = 0; j < a.cols(); ++j) {
      for (Eigen::Index i = 0; i < a.rows(); ++i) {
        a.adj().coeffRef(i, j) -= vi.adj().coeff(i, j);
      }
    }
  });
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_REV_FUN_FABS_HPP
#define STAN_MATH_REV_FUN_FABS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Return the absolute value of the variable (cmath).
 *
 * Choosing an arbitrary value at the non-differentiable point 0,
 *
 * \f$\frac{d}{dx}|x| = \mbox{sgn}(x)\f$.
 *
 * where \f$\mbox{sgn}(x)\f$ is the signum function, taking values
 * -1 if \f$x < 0\f$, 0 if \f$x == 0\f$, and 1 if \f$x == 1\f$.
 *
 * The function <code>abs()</code> provides the same behavior, with
 * <code>abs()</code> defined in stdlib.h and <code>fabs()</code>
 * defined in <code>cmath</code>.
 * The derivative is 0 if the input is 0.
 *
 * Returns std::numeric_limits<double>::quiet_NaN() for NaN inputs.
 *
 *
   \f[
   \mbox{fabs}(x) =
   \begin{cases}
     |x| & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fabs}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } x < 0 \\
     0 & \mbox{if } x = 0 \\
     1 & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Input variable.
 * @return Absolute value of variable.
 */
inline var fabs(const var& a) {
  if (a.val() > 0.0) {
    return a;
  } else if (a.val() < 0.0) {
    return -a;
  } else if (a.val() == 0) {
    return var(new vari(0));
  } else {
    return make_callback_var(NOT_A_NUMBER,
                             [a](auto& vi) mutable { a.adj() = NOT_A_NUMBER; });
  }
}

/**
 * Return the absolute value of the variable (cmath).
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param a Input variable.
 * @return Absolute value of variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto fabs(const VarMat& a) {
  return make_callback_var(
      a.val().unaryExpr([](const auto x) {
        if (unlikely(is_nan(x))) {
          return NOT_A_NUMBER;
        } else {
          return std::abs(x);
        }
      }),
      [a](const auto& vi) mutable {
        for (Eigen::Index j = 0; j < vi.cols(); ++j) {
          for (Eigen::Index i = 0; i < vi.rows(); ++i) {
            const auto x = a.val().coeffRef(i, j);
            if (unlikely(is_nan(x))) {
              a.adj().coeffRef(i, j) = NOT_A_NUMBER;
            } else if (x < 0.0) {
              a.adj().coeffRef(i, j) -= vi.adj_.coeff(i, j);
            } else {
              a.adj().coeffRef(i, j) += vi.adj_.coeff(i, j);
            }
          }
        }
      });
}

}  // namespace math
}  // namespace stan
#endif

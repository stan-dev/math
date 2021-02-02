#ifndef STAN_MATH_OPENCL_REV_fmod_HPP
#define STAN_MATH_OPENCL_REV_fmod_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Return the floating point remainder after dividing the
 * first variable by the second (cmath).
 *
 * The partial derivatives with respect to the variables are defined
 * everywhere but where \f$x = y\f$, but we set these to match other values,
 * with
 *
 * \f$\frac{\partial}{\partial x} \mbox{fmod}(x, y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} \mbox{fmod}(x, y) = -\lfloor \frac{x}{y}
 \rfloor\f$.
 *
 *
   \f[
   \mbox{fmod}(x, y) =
   \begin{cases}
     x - \lfloor \frac{x}{y}\rfloor y & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fmod}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fmod}(x, y)}{\partial y} =
   \begin{cases}
     -\lfloor \frac{x}{y}\rfloor & \mbox{if } -\infty\leq x, y\leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable.
 * @param b Second variable.
 * @return Floating pointer remainder of dividing the first variable
 * by the second.
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> fmod(T_a&& a, T_b&& b) {
  using std::isnan;
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = fmod(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto any_nan = isnan(value_of(a_arena)) || isnan(value_of(b_arena));
        auto a_is_max = value_of(a_arena) > value_of(b_arena);
        auto a_deriv = select(any_nan, NOT_A_NUMBER, res.adj());
        auto b_deriv = select(
            any_nan, NOT_A_NUMBER,
            elt_multiply(-res.adj(), trunc(elt_divide(value_of(a_arena),
                                                      value_of(b_arena)))));

        adjoint_results(a_arena, b_arena) += expressions(a_deriv, b_deriv);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

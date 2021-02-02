#ifndef STAN_MATH_OPENCL_REV_FDIM_HPP
#define STAN_MATH_OPENCL_REV_FDIM_HPP
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
 * Return the positive difference between the first variable's the value
 * and the second's (C99, C++11).
 *
 * The function values and derivatives are defined by
 *
   \f[
   \mbox{fdim}(x, y) =
   \begin{cases}
     x-y & \mbox{if } x  > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fdim}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } x > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fdim}(x, y)}{\partial y} =
   \begin{cases}
    -1 & \mbox{if } x > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable.
 * @param b Second variable.
 * @return The positive difference between the first and second
 * variable.
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> fdim(T_a&& a, T_b&& b) {
  using std::isnan;
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = fdim(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto nan_check
            = select(isnan(value_of(a_arena)) || isnan(value_of(b_arena)),
                     NOT_A_NUMBER, 0.0);
        auto a_is_max = value_of(a_arena) > value_of(b_arena);
        auto a_deriv = select(a_is_max, res.adj(), nan_check);
        auto b_deriv = select(a_is_max, -res.adj(), nan_check);

        adjoint_results(a_arena, b_arena) += expressions(a_deriv, b_deriv);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

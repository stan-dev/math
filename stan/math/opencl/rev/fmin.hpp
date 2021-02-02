#ifndef STAN_MATH_OPENCL_REV_FMIN_HPP
#define STAN_MATH_OPENCL_REV_FMIN_HPP
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
 * Return the lesser of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @param a First argument.
 * @param b Second argument.
 * @return Minimum of a or b and if one is NaN return the other
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> fmin(T_a&& a, T_b&& b) {
  using std::isnan;
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = fmin(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto both_nan = isnan(value_of(a_arena)) && isnan(value_of(b_arena));
        auto a_is_min = value_of(a_arena) < value_of(b_arena);
        auto a_deriv
            = select(both_nan, NOT_A_NUMBER, select(a_is_min, res.adj(), 0.0));
        auto b_deriv
            = select(both_nan, NOT_A_NUMBER, select(!a_is_min, res.adj(), 0.0));

        adjoint_results(a_arena, b_arena) += expressions(a_deriv, b_deriv);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

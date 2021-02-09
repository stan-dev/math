#ifndef STAN_MATH_OPENCL_REV_LOG_INV_LOGIT_DIFF_HPP
#define STAN_MATH_OPENCL_REV_LOG_INV_LOGIT_DIFF_HPP
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
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Returns the natural logarithm of the difference of the
 * inverse logits of the specified arguments.
 *
 * @tparam T_x type of x argument
 * @tparam T_y type of y argument
 * @param x first argument
 * @param y second argument
 * @return Result of log difference of inverse logits of arguments.
 */
template <typename T_x, typename T_y,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_y>* = nullptr,
          require_any_var_t<T_x, T_y>* = nullptr,
          require_any_not_stan_scalar_t<T_x, T_y>* = nullptr>
inline var_value<matrix_cl<double>> log_inv_logit_diff(T_x&& x, T_y&& y) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_y> y_arena = std::forward<T_y>(y);

  return make_callback_var(
      log_inv_logit_diff(value_of(x_arena), value_of(y_arena)),
      [x_arena, y_arena](const vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(x_arena, y_arena) += expressions(
            -elt_multiply(res.adj(),
                          inv(expm1(value_of(y_arena) - value_of(x_arena)))
                              + inv_logit(value_of(x_arena))),
            -elt_multiply(res.adj(),
                          inv(expm1(value_of(x_arena) - value_of(y_arena)))
                              + inv_logit(value_of(y_arena))));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

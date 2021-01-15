#ifndef STAN_MATH_OPENCL_REV_LOG_DIFF_EXP_HPP
#define STAN_MATH_OPENCL_REV_LOG_DIFF_EXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns  the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 *
 * @tparam T_x type of x argument
 * @tparam T_y type of y argument
 * @param x first argument
 * @param y second argument
 * @return Result the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 */
template <
    typename T_x, typename T_y,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_x, T_y>* = nullptr,
    require_any_var_t<T_x, T_y>* = nullptr>
inline var_value<matrix_cl<double>> log_diff_exp(T_x&& x, T_y&& y) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_y> y_arena = std::forward<T_y>(y);

  matrix_cl<double> res_val
      = log_diff_exp(value_of(x_arena), value_of(y_arena));

  return make_callback_var(
      res_val,
      [x_arena, y_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_x>::value && !is_constant<T_y>::value) {
          auto x_deriv = -elt_divide(
              res.adj(), expm1(value_of(y_arena) - value_of(x_arena)));
          auto y_deriv = -elt_divide(
              res.adj(), expm1(value_of(x_arena) - value_of(y_arena)));

          results(adjoint_of(x_arena), adjoint_of(y_arena))
              += expressions(calc_if<is_var<T_x>::value>(x_deriv),
                             calc_if<is_var<T_y>::value>(y_deriv));
        } else if (!is_constant<T_x>::value) {
          auto& x_adj = forward_as<var_value<matrix_cl<double>>>(x_arena).adj();
          x_adj = x_adj
                  + select(res.val() == NEGATIVE_INFTY,
                           select(value_of(y_arena) == NEGATIVE_INFTY,
                                  res.adj(), INFTY * res.adj()),
                           -elt_divide(res.adj(), expm1(value_of(y_arena)
                                                        - value_of(x_arena))));
        } else {
          auto& y_adj = forward_as<var_value<matrix_cl<double>>>(y_arena).adj();
          y_adj = y_adj
                  + select(res.val() == NEGATIVE_INFTY, -(res.adj() * INFTY),
                           -elt_divide(res.adj(), expm1(value_of(y_arena)
                                                        - value_of(x_arena))));
        }
      });
}

/**
 * Returns  the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 *
 * @tparam T_x type of kernel generator expression for x
 * @tparam T_y type of scalar y argument
 * @param x kernel generator expression
 * @param y scalar
 * @return Result the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 */
template <typename T_x, typename T_y,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr,
          require_stan_scalar_t<T_y>* = nullptr,
          require_any_var_t<T_x, T_y>* = nullptr>
inline var_value<matrix_cl<double>> log_diff_exp(T_x&& x, const T_y& y) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);

  matrix_cl<double> res_val = log_diff_exp(value_of(x_arena), value_of(y));

  return make_callback_var(
      res_val, [x_arena, y](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_x>::value) {
          auto& x_adj = forward_as<var_value<matrix_cl<double>>>(x_arena).adj();
          x_adj
              = x_adj
                - elt_divide(res.adj(), expm1(value_of(y) - value_of(x_arena)));
        }
        if (!is_constant<T_y>::value) {
          adjoint_of(y) -= sum(
              elt_divide(res.adj(), expm1(value_of(x_arena) - value_of(y))));
        }
      });
}

/**
 * Returns the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 *
 * @tparam T_x type of scalar x argument
 * @tparam T_y type of kernel generator expression for y
 * @param x kernel generator expression
 * @param y scalar
 * @return Result the natural logarithm of the difference
 * of the natural exponentiation of x and
 * the natural exponentiation of y.
 */
template <typename T_x, typename T_y,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr,
          require_stan_scalar_t<T_x>* = nullptr,
          require_any_var_t<T_x, T_y>* = nullptr>
inline var_value<matrix_cl<double>> log_diff_exp(const T_x& x, T_y&& y) {
  arena_t<T_y> y_arena = std::forward<T_y>(y);

  matrix_cl<double> res_val = log_diff_exp(value_of(x), value_of(y_arena));

  return make_callback_var(
      res_val, [x, y_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_x>::value) {
          adjoint_of(x) -= sum(
              elt_divide(res.adj(), expm1(value_of(y_arena) - value_of(x))));
        }
        if (!is_constant<T_y>::value) {
          auto& y_adj = forward_as<var_value<matrix_cl<double>>>(y_arena).adj();
          y_adj
              = y_adj
                - elt_divide(res.adj(), expm1(value_of(x) - value_of(y_arena)));
        }
      });
}
}  // namespace math
}  // namespace stan

#endif
#endif

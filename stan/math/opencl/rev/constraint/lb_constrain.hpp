#ifndef STAN_MATH_OPENCL_REV_CONSTRAINT_LB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_REV_CONSTRAINT_LB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T_x type of unconstrained input
 * @tparam T_lb type of lower bound
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @return constrained matrix
 */
template <typename T_x, typename T_lb,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_lb>* = nullptr,
          require_any_var_t<T_x, T_lb>* = nullptr,
          require_any_not_stan_scalar_t<T_x, T_lb>* = nullptr>
inline var_value<matrix_cl<double>> lb_constrain(T_x&& x, T_lb&& lb) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_lb> lb_arena = std::forward<T_lb>(lb);

  return make_callback_var(
      lb_constrain(value_of(x_arena), value_of(lb_arena)),
      [x_arena, lb_arena](vari_value<matrix_cl<double>>& res) mutable {
        auto lb_inf = value_of(lb_arena) == NEGATIVE_INFTY;
        adjoint_results(x_arena, lb_arena) += expressions(
            select(lb_inf, res.adj(),
                   elt_multiply(res.adj(), exp(value_of(x_arena)))),
            select(lb_inf, 0, res.adj()));
      });
}

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T_x type of unconstrained input
 * @tparam T_lb type of lower bound
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @param[in,out] lp reference to log probability to increment
 * @return constrained matrix
 */
template <typename T_x, typename T_lb,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_lb>* = nullptr,
          require_any_var_t<T_x, T_lb>* = nullptr,
          require_any_not_stan_scalar_t<T_x, T_lb>* = nullptr>
inline var_value<matrix_cl<double>> lb_constrain(T_x&& x, T_lb&& lb, var& lp) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_lb> lb_arena = std::forward<T_lb>(lb);

  double lp_inc = 0;
  matrix_cl<double> res
      = lb_constrain(value_of(x_arena), value_of(lb_arena), lp_inc);
  lp += lp_inc;

  return make_callback_var(
      std::move(res),
      [x_arena, lb_arena, lp](vari_value<matrix_cl<double>>& res) mutable {
        auto lb_inf = value_of(lb_arena) == NEGATIVE_INFTY;
        adjoint_results(x_arena, lb_arena) += expressions(
            select(lb_inf, res.adj(),
                   elt_multiply(res.adj(), exp(value_of(x_arena))) + lp.adj()),
            select(lb_inf, 0, res.adj()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

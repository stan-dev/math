#ifndef STAN_MATH_OPENCL_REV_CONSTRAINT_LUB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_REV_CONSTRAINT_LUB_CONSTRAIN_HPP
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
 * Return the bounded value for the specified unconstrained input
 * and specified bounds.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T_x type of the unconstrained input
 * @tparam T_lb type of lower bound
 * @tparam T_ub type of upper bound
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @param[in] ub upper bound
 * @return constrained matrix
 */
template <
    typename T_x, typename T_lb, typename T_ub,
    require_all_prim_or_rev_kernel_expression_t<T_x, T_lb, T_ub>* = nullptr,
    require_any_var_t<T_x, T_lb, T_ub>* = nullptr,
    require_any_not_stan_scalar_t<T_x, T_lb, T_ub>* = nullptr>
inline var_value<matrix_cl<double>> lub_constrain(T_x&& x, T_lb&& lb,
                                                  T_ub&& ub) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_lb> lb_arena = std::forward<T_lb>(lb);
  arena_t<T_ub> ub_arena = std::forward<T_ub>(ub);

  return make_callback_var(
      lub_constrain(value_of(x_arena), value_of(lb_arena), value_of(ub_arena)),
      [x_arena, lb_arena,
       ub_arena](vari_value<matrix_cl<double>>& res) mutable {
        auto lb_inf = value_of(lb_arena) == NEGATIVE_INFTY;
        auto ub_inf = value_of(ub_arena) == INFTY;
        auto inv_logit_x = inv_logit(value_of(x_arena));
        auto one_m_inv_logit_x = 1.0 - inv_logit_x;
        auto exp_x = exp(value_of(x_arena));
        auto res_adj_exp_x = elt_multiply(res.adj(), exp_x);
        adjoint_results(x_arena, lb_arena, ub_arena) += expressions(
            select(lb_inf, select(ub_inf, res.adj(), -res_adj_exp_x),
                   select(ub_inf, res_adj_exp_x,
                          elt_multiply(
                              elt_multiply(res.adj(), (value_of(ub_arena)
                                                       - value_of(lb_arena))),
                              elt_multiply(inv_logit_x, one_m_inv_logit_x)))),
            select(lb_inf, 0,
                   select(ub_inf, res.adj(),
                          elt_multiply(res.adj(), one_m_inv_logit_x))),
            select(ub_inf, 0,
                   select(lb_inf, res.adj(),
                          elt_multiply(res.adj(), inv_logit_x))));
      });
}

/**
 * Return the bounded value for the specified unconstrained input
 * and specified bounds.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T_x type of the unconstrained input
 * @tparam T_lb type of lower bound
 * @tparam T_ub type of upper bound
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @param[in] ub upper bound
 * @param[in,out] lp reference to log probability to increment
 * @return constrained matrix
 */
template <
    typename T_x, typename T_lb, typename T_ub,
    require_all_prim_or_rev_kernel_expression_t<T_x, T_lb, T_ub>* = nullptr,
    require_any_var_t<T_x, T_lb, T_ub>* = nullptr,
    require_any_not_stan_scalar_t<T_x, T_lb, T_ub>* = nullptr>
inline var_value<matrix_cl<double>> lub_constrain(T_x&& x, T_lb&& lb, T_ub&& ub,
                                                  var& lp) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_lb> lb_arena = std::forward<T_lb>(lb);
  arena_t<T_ub> ub_arena = std::forward<T_ub>(ub);

  double lp_inc = 0;
  matrix_cl<double> res = lub_constrain(value_of(x_arena), value_of(lb_arena),
                                        value_of(ub_arena), lp_inc);
  lp += lp_inc;

  return make_callback_var(
      std::move(res), [x_arena, lb_arena, ub_arena,
                       lp](vari_value<matrix_cl<double>>& res) mutable {
        auto lb_inf = value_of(lb_arena) == NEGATIVE_INFTY;
        auto ub_inf = value_of(ub_arena) == INFTY;
        auto inv_logit_x = inv_logit(value_of(x_arena));
        auto one_m_inv_logit_x = 1.0 - inv_logit_x;
        auto diff = value_of(ub_arena) - value_of(lb_arena);
        auto one_over_diff = elt_divide(1.0, diff);
        auto exp_x = exp(value_of(x_arena));
        auto res_adj_exp_x = elt_multiply(res.adj(), exp_x);
        adjoint_results(x_arena, lb_arena, ub_arena) += expressions(
            select(lb_inf, select(ub_inf, res.adj(), lp.adj() - res_adj_exp_x),
                   select(ub_inf, res_adj_exp_x + lp.adj(),
                          elt_multiply(
                              elt_multiply(res.adj(), diff),
                              elt_multiply(inv_logit_x, one_m_inv_logit_x))
                              + lp.adj() * (1.0 - 2.0 * inv_logit_x))),
            select(lb_inf, 0.0,
                   select(ub_inf, res.adj(),
                          elt_multiply(res.adj(), one_m_inv_logit_x)
                              - one_over_diff * lp.adj())),
            select(ub_inf, 0.0,
                   select(lb_inf, res.adj(),
                          elt_multiply(res.adj(), inv_logit_x)
                              + one_over_diff * lp.adj())));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

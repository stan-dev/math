#ifndef STAN_MATH_OPENCL_REV_CONSTRAINT_UB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_REV_CONSTRAINT_UB_CONSTRAIN_HPP
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
 * Return the upper-bounded value for the specified unconstrained input
 * and specified upper bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the constant upper bound.
 *
 * @tparam T_x kernel generator expression
 * @tparam T_ub kernel generator expression
 * @param[in] x unconstrained input
 * @param[in] ub upper bound
 * @return constrained matrix
 */
template <typename T_x, typename T_ub,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_ub>* = nullptr,
          require_any_var_t<T_x, T_ub>* = nullptr,
          require_any_not_stan_scalar_t<T_x, T_ub>* = nullptr>
inline var_value<matrix_cl<double>> ub_constrain(T_x&& x, T_ub&& ub) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_ub> ub_arena = std::forward<T_ub>(ub);

  return make_callback_var(
      ub_constrain(value_of(x_arena), value_of(ub_arena)),
      [x_arena, ub_arena](vari_value<matrix_cl<double>>& res) mutable {
        auto ub_inf = value_of(ub_arena) == INFTY;
        adjoint_results(x_arena, ub_arena) += expressions(
            select(ub_inf, res.adj(),
                   -elt_multiply(res.adj(), exp(value_of(x_arena)))),
            select(ub_inf, 0, res.adj()));
      });
}

/**
 * Return the upper-bounded value for the specified unconstrained input
 * and specified upper bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the constant upper bound.
 *
 * @tparam T_x kernel generator expression
 * @tparam T_ub kernel generator expression
 * @param[in] x unconstrained input
 * @param[in] ub upper bound
 * @param[in,out] lp reference to log probability to increment
 * @return constrained matrix
 */
template <typename T_x, typename T_ub,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_ub>* = nullptr,
          require_any_var_t<T_x, T_ub>* = nullptr,
          require_any_not_stan_scalar_t<T_x, T_ub>* = nullptr>
inline var_value<matrix_cl<double>> ub_constrain(T_x&& x, T_ub&& ub, var& lp) {
  arena_t<T_x> x_arena = std::forward<T_x>(x);
  arena_t<T_ub> ub_arena = std::forward<T_ub>(ub);

  double lp_inc = 0;
  matrix_cl<double> res
      = ub_constrain(value_of(x_arena), value_of(ub_arena), lp_inc);
  lp += lp_inc;

  return make_callback_var(
      std::move(res),
      [x_arena, ub_arena, lp](vari_value<matrix_cl<double>>& res) mutable {
        auto ub_inf = value_of(ub_arena) == INFTY;
        adjoint_results(x_arena, ub_arena) += expressions(
            select(ub_inf, res.adj(),
                   -elt_multiply(res.adj(), exp(value_of(x_arena))) + lp.adj()),
            select(ub_inf, 0, res.adj()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

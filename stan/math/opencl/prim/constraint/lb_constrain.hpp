#ifndef STAN_MATH_OPENCL_PRIM_CONSTRAINT_LB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_CONSTRAINT_LB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

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
 * @tparam T kernel generator expression
 * @tparam L kernel generator expression
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @return constrained matrix
 */
template <typename T, typename L,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<L>* = nullptr>
inline auto lb_constrain(T&& x, L&& lb) {
  return make_holder_cl(
      [](auto& x_, auto& lb_) {
        return select(lb_ == NEGATIVE_INFTY, x_, lb_ + exp(x_));
      },
      std::forward<T>(x), std::forward<L>(lb));
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
 * @tparam T kernel generator expression
 * @tparam L kernel generator expression
 * @param[in] x unconstrained input
 * @param[in] lb lower bound
 * @param[in,out] lp reference to log probability to increment
 * @return constrained matrix
 */
template <typename T, typename L,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  matrix_cl<double> lp_inc;
  matrix_cl<double> res;
  auto lb_inf = lb == NEGATIVE_INFTY;
  auto lp_inc_expr = sum_2d(select(lb_inf, 0.0, x));
  auto res_expr = select(lb_inf, x, lb + exp(x));
  results(lp_inc, res) = expressions(lp_inc_expr, res_expr);
  lp += sum(from_matrix_cl(lp_inc));
  return res;
}

}  // namespace math
}  // namespace stan
#endif
#endif

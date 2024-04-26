#ifndef STAN_MATH_OPENCL_PRIM_CONSTRAINT_LUB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_CONSTRAINT_LUB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Return the lower and upper-bounded matrix derived by
 * transforming the specified free matrix given the specified
 * lower and upper bounds.
 *
 * <p>The transform is the transformed and scaled inverse logit,
 *
 * <p>\f$f(x) = L + (U - L) \mbox{logit}^{-1}(x)\f$
 *
 * @tparam T matrix expression type
 * @tparam L lower bound expression type
 * @tparam U upper bound expression type
 * @param[in] x Free matrix to transform.
 * @param[in] lb Lower bound
 * @param[in] ub Upper bound
 * @return Lower- and upper-bounded matrix derived from transforming
 *   the free matrix.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<L, U>* = nullptr>
inline matrix_cl<double> lub_constrain(const T& x, const L& lb, const U& ub) {
  auto diff = ub - lb;
  auto lb_inf = lb == NEGATIVE_INFTY;
  auto ub_inf = ub == INFTY;

  auto check
      = check_cl("lub_constrain (OpenCL)", "(ub - lb)", diff, "positive");
  matrix_cl<double> res;

  results(check, res) = expressions(
      diff > 0.0, select(lb_inf, select(ub_inf, x, ub - exp(x)),
                         select(ub_inf, exp(x) + lb,
                                elt_multiply(diff, inv_logit(x)) + lb)));
  return res;
}

/**
 * Return the lower and upper-bounded matrix derived by
 * transforming the specified free matrix given the specified
 * lower and upper bounds.
 *
 * <p>The transform is the transformed and scaled inverse logit,
 *
 * <p>\f$f(x) = L + (U - L) \mbox{logit}^{-1}(x)\f$
 *
 * @tparam T matrix expression type
 * @tparam L lower bound expression type
 * @tparam U upper bound expression type
 * @param[in] x Free matrix to transform.
 * @param[in] lb Lower bound
 * @param[in] ub Upper bound
 * @param[in,out] lp Log probability scalar reference
 * @return Lower- and upper-bounded matrix derived from transforming
 *   the free matrix.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<L, U>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  auto diff = ub - lb;
  auto lb_inf = lb == NEGATIVE_INFTY;
  auto ub_inf = ub == INFTY;
  auto abs_x = fabs(x);
  auto check
      = check_cl("lub_constrain (OpenCL)", "(ub - lb)", diff, "positive");

  matrix_cl<double> res;
  matrix_cl<double> lp_inc;

  auto lp_inc_expr = sum_2d(
      select(lb_inf, select(ub_inf, 0.0, x),
             select(ub_inf, x, log(diff) - abs_x - 2.0 * log1p_exp(-abs_x))));
  auto res_expr = select(
      lb_inf, select(ub_inf, x, ub - exp(x)),
      select(ub_inf, exp(x) + lb, elt_multiply(diff, inv_logit(x)) + lb));

  results(check, res, lp_inc) = expressions(diff > 0.0, res_expr, lp_inc_expr);

  lp += sum(from_matrix_cl(lp_inc));

  return res;
}

}  // namespace math
}  // namespace stan
#endif
#endif

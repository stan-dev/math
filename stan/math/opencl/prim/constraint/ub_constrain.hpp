#ifndef STAN_MATH_OPENCL_PRIM_CONSTRAINT_UB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_CONSTRAINT_UB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/opencl/scalar_cl_reduce.hpp>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<U>* = nullptr>
inline auto ub_constrain(T&& x, U&& ub) {
  return make_holder_cl(
      [](auto&& x_, auto&& ub_) {
        // the operands are moved into the returned expression, so device
        // scalars must not be referenced through locals of this lambda
        return select(opencl::internal::as_operand(ub_) == INFTY, x_,
                      opencl::internal::as_operand(ub_) - exp(x_));
      },
      std::forward<T>(x), std::forward<U>(ub));
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @param[in,out] lp reference to log probability to increment
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U, typename T_lp,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<U>* = nullptr,
          require_t<math::disjunction<std::is_same<T_lp, double>,
                                      is_prim_scalar_cl<T_lp>>>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub_in, T_lp& lp) {
  const auto& ub = opencl::internal::as_operand(ub_in);
  matrix_cl<double> lp_inc;
  matrix_cl<double> res;
  auto ub_inf = ub == INFTY;
  auto lp_inc_expr = sum_2d(select(ub_inf, 0.0, x));
  auto res_expr = select(ub_inf, x, ub - exp(x));
  results(lp_inc, res) = expressions(lp_inc_expr, res_expr);
  opencl::internal::add_sum(lp, lp_inc);
  return res;
}

}  // namespace math
}  // namespace stan
#endif
#endif

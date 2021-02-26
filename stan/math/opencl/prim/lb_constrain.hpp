#ifndef STAN_MATH_OPENCL_PRIM_LB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_LB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

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
template <typename T, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<U>* = nullptr>
inline auto lb_constrain(T&& x, U&& lb) {
  return make_holder_cl(
      [](auto& x_, auto& lb_) {
        return select(lb_ == NEGATIVE_INFTY, x_, lb_ + exp(x_));
      },
      std::forward<T>(x), std::forward<U>(lb));
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
template <typename T, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<U>* = nullptr>
inline auto lb_constrain(T&& x, U&& lb, return_type_t<T, U>& lp) {
  lp += sum(x);
  return lb_constrain(std::forward<T>(x), std::forward<U>(lb));
}

}  // namespace math
}  // namespace stan
#endif
#endif

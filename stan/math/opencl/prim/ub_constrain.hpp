#ifndef STAN_MATH_OPENCL_PRIM_UB_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_UB_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

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
      [](auto& x_, auto& ub_) {
        return select(ub_ == INFTY, x_, ub_ - exp(x_));
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
template <typename T, typename U,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr,
          require_all_kernel_expressions_t<U>* = nullptr>
inline auto ub_constrain(T&& x, U&& ub, return_type_t<T, U>& lp) {
  lp += sum(x);
  return ub_constrain(std::forward<T>(x), std::forward<U>(ub));
}

}  // namespace math
}  // namespace stan
#endif
#endif

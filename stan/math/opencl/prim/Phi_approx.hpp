#ifndef STAN_MATH_OPENCL_PRIM_PHI_APPROX_HPP
#define STAN_MATH_OPENCL_PRIM_PHI_APPROX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise logistic approximation to the standard normal CDF,
 * `inv_logit(PHI_APPROX_CUBIC_COEF * x^3 + PHI_APPROX_LINEAR_COEF * x)`, of a
 * kernel generator expression.
 *
 * @tparam T type of the kernel generator expression
 * @param x argument
 * @return `Phi_approx` of each element of x
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto Phi_approx(T&& x) {
  return make_holder_cl(
      [](auto&& a) {
        return inv_logit(PHI_APPROX_CUBIC_COEF * pow(a, 3.0)
                         + PHI_APPROX_LINEAR_COEF * a);
      },
      to_ref(std::forward<T>(x)));
}

}  // namespace math
}  // namespace stan

#endif
#endif

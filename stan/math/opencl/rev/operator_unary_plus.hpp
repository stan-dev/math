#ifndef STAN_MATH_OPENCL_REV_OPERATOR_UNARY_PLUS_HPP
#define STAN_MATH_OPENCL_REV_OPERATOR_UNARY_PLUS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the unary plus of the input.
 *
 * @param M input kernel expression
 * @return result of unary plus of the input.
 */
template <typename T,
          require_var_vt<is_kernel_expression_and_not_scalar, T>* = nullptr>
inline T operator+(T&& M) {
  return std::forward<T>(M);
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_TRUNC_HPP
#define STAN_MATH_OPENCL_REV_TRUNC_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `trunc()` of the
 * input `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `trunc()` of the input argument.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> trunc(const var_value<T>& A) {
  return var_value<matrix_cl<double>>(trunc(A.val()));
}

}  // namespace math
}  // namespace stan

#endif
#endif

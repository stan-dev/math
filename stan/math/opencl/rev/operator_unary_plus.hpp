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
inline var_value<matrix_cl<double>> operator+(
    const var_value<matrix_cl<double>>& M) {
  return M;
}

}  // namespace math
}  // namespace stan

#endif
#endif

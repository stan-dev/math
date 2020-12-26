#ifndef STAN_MATH_OPENCL_REV_OPERATOR_UNARY_MINUS_HPP
#define STAN_MATH_OPENCL_REV_OPERATOR_UNARY_MINUS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the unary minus of the input.
 *
 * @param M input kernel expression
 * @return result of unary minus of the input.
 */
inline var_value<matrix_cl<double>> operator-(
    const var_value<matrix_cl<double>>& M) {
  var_value<matrix_cl<double>> res = -M.val();

  reverse_pass_callback([M, res]() mutable { M.adj() = M.adj() - res.adj(); });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_EXP2_HPP
#define STAN_MATH_OPENCL_REV_EXP2_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `exp2()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `exp2()` of the input.
 */
inline var_value<matrix_cl<double>> exp2(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = exp2(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj() + elt_multiply(elt_multiply(res.adj(), res.val()), LOG_TWO);
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_ACOS_HPP
#define STAN_MATH_OPENCL_REV_ACOS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `acos()` of a var_value<matrix_cl<double>>
 * in radians.
 *
 * @param A argument
 * @return Elementwise `acos()` of the input, in radians.
 */
inline var_value<matrix_cl<double>> acos(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = acos(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          - elt_divide(res.adj(), sqrt(1.0 - elt_multiply(A.val(), A.val())));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_ASINH_HPP
#define STAN_MATH_OPENCL_REV_ASINH_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `asinh()` of a var_value<matrix_cl<double>>
 * in radians.
 *
 * @param A argument
 * @return Elementwise `asinh()` of the input, in radians.
 */
inline var_value<matrix_cl<double>> asinh(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = asinh(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_divide(res.adj(), sqrt(elt_multiply(A.val(), A.val()) + 1.0));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

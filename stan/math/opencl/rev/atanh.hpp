#ifndef STAN_MATH_OPENCL_REV_ATANH_HPP
#define STAN_MATH_OPENCL_REV_ATANH_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `atanh()` of a var_value<matrix_cl<double>>
 * in radians.
 *
 * @param A argument
 * @return Elementwise `atanh()` of the input, in radians.
 */
inline var_value<matrix_cl<double>> atanh(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = atanh(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + elt_divide(res.adj(), (1.0 - elt_multiply(A.val(), A.val())));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

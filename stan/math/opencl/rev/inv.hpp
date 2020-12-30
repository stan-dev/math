#ifndef STAN_MATH_OPENCL_REV_INV_HPP
#define STAN_MATH_OPENCL_REV_INV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `inv()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `inv()` of the input.
 */
inline var_value<matrix_cl<double>> inv(const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = inv(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() - elt_divide(res.adj(), elt_multiply(A.val(), A.val()));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

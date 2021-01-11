#ifndef STAN_MATH_OPENCL_REV_SQUARE_HPP
#define STAN_MATH_OPENCL_REV_SQUARE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `square()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `square()` of the input.
 */
inline var_value<matrix_cl<double>> square(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = square(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj() = A.adj() + 2.0 * elt_multiply(res.adj(), A.val());
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

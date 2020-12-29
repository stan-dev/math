#ifndef STAN_MATH_OPENCL_REV_ERF_HPP
#define STAN_MATH_OPENCL_REV_ERF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `erf()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `erf()` of the input.
 */
inline var_value<matrix_cl<double>> erf(const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = erf(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_multiply(res.adj(),
                         elt_multiply(TWO_OVER_SQRT_PI,
                                      exp(-elt_multiply(A.val(), A.val()))));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

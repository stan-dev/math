#ifndef STAN_MATH_OPENCL_REV_CBRT_HPP
#define STAN_MATH_OPENCL_REV_CBRT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `cbrt()` of the input
 * `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `cbrt()` of the input argument.
 */
inline var_value<matrix_cl<double>> cbrt(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = cbrt(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_divide(res.adj(),
                       elt_multiply(3.0, elt_multiply(res.val(), res.val())));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

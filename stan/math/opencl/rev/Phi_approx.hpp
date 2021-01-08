#ifndef STAN_MATH_OPENCL_REV_PHI_APPROX_HPP
#define STAN_MATH_OPENCL_REV_PHI_APPROX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `Phi_approx()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `Phi_approx()` of the input.
 */
inline var_value<matrix_cl<double>> Phi_approx(
    const var_value<matrix_cl<double>>& A) {
  var_value<matrix_cl<double>> res = Phi_approx(A.val());

  reverse_pass_callback([A, res]() mutable {
    A.adj()
        = A.adj()
          + elt_multiply(
                res.adj(),
                elt_multiply(
                    res.val(),
                    elt_multiply((1 - res.val()),
                                 (3.0 * 0.07056 * elt_multiply(A.val(), A.val())
                                  + 1.5976))));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

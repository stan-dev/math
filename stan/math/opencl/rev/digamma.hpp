#ifndef STAN_MATH_OPENCL_REV_DIGAMMA_HPP
#define STAN_MATH_OPENCL_REV_DIGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `digamma()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `digamma()` of the input.
 */
inline var_value<matrix_cl<double>> digamma(
    const var_value<matrix_cl<double>>& A) {
  return make_callback_var(
      digamma(A.val()), [A](const vari_value<matrix_cl<double>> res) mutable {
        A.adj() += elt_multiply(res.adj(), trigamma(A.val()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

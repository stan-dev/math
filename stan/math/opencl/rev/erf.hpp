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
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> erf(const var_value<T>& A) {
  return make_callback_var(
      erf(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(
            res.adj(), elt_multiply(TWO_OVER_SQRT_PI, exp(-square(A.val()))));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

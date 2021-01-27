#ifndef STAN_MATH_OPENCL_REV_LGAMMA_HPP
#define STAN_MATH_OPENCL_REV_LGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `lgamma()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `lgamma()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> lgamma(const var_value<T>& A) {
  return make_callback_var(
      lgamma(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(res.adj(), digamma(A.val()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

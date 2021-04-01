#ifndef STAN_MATH_OPENCL_REV_INV_CLOGLOG_HPP
#define STAN_MATH_OPENCL_REV_INV_CLOGLOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `inv_cloglog()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `inv_cloglog()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> inv_cloglog(const var_value<T>& A) {
  return make_callback_var(
      inv_cloglog(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(res.adj(), exp(A.val() - exp(A.val())));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

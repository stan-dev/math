#ifndef STAN_MATH_OPENCL_REV_LOGIT_HPP
#define STAN_MATH_OPENCL_REV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `logit()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `logit()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> logit(const var_value<T>& A) {
  return make_callback_var(
      logit(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_divide(res.adj(), A.val() - square(A.val()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

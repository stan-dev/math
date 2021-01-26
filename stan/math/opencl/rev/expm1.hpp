#ifndef STAN_MATH_OPENCL_REV_EXPM1_HPP
#define STAN_MATH_OPENCL_REV_EXPM1_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `expm1()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `expm1()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> expm1(const var_value<T>& A) {
  return make_callback_var(
      expm1(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(res.adj(), res.val() + 1.0);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

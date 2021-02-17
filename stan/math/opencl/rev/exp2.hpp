#ifndef STAN_MATH_OPENCL_REV_EXP2_HPP
#define STAN_MATH_OPENCL_REV_EXP2_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `exp2()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `exp2()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> exp2(const var_value<T>& A) {
  return make_callback_var(
      exp2(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(elt_multiply(res.adj(), res.val()), LOG_TWO);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

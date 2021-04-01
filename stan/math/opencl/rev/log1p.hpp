#ifndef STAN_MATH_OPENCL_REV_LOG1P_HPP
#define STAN_MATH_OPENCL_REV_LOG1P_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `log1p()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `log1p()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> log1p(const var_value<T>& A) {
  return make_callback_var(log1p(A.val()),
                           [A](vari_value<matrix_cl<double>>& res) mutable {
                             A.adj() += elt_divide(res.adj(), 1.0 + A.val());
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_CEIL_HPP
#define STAN_MATH_OPENCL_REV_CEIL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `ceil()` of the specified variable.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `ceil()` of the input argument.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> ceil(const var_value<T>& A) {
  return make_callback_var(
      ceil(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() = select(isnan(A.val()),
                         constant(NOT_A_NUMBER, A.rows(), A.cols()), A.adj());
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

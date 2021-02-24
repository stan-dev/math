#ifndef STAN_MATH_OPENCL_REV_FABS_HPP
#define STAN_MATH_OPENCL_REV_FABS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `fabs()` of the
 * input `var_value<matrix_cl<double>>`.
 *
 * @param A input `var_value<matrix_cl<double>>`.
 * @return Elementwise `fabs()` of the input argument.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> fabs(const var_value<T>& A) {
  return make_callback_var(
      fabs(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() = select(
            isnan(A.val()), NOT_A_NUMBER,
            select(A.val() < 0.0, A.adj() - res.adj(), A.adj() + res.adj()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

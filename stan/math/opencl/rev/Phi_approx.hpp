#ifndef STAN_MATH_OPENCL_REV_PHI_APPROX_HPP
#define STAN_MATH_OPENCL_REV_PHI_APPROX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `Phi_approx()` of a var_value<matrix_cl<double>>.
 *
 * @param A argument
 * @return Elementwise `Phi_approx()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> Phi_approx(const var_value<T>& A) {
  return make_callback_var(
      Phi_approx(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(
            elt_multiply(elt_multiply(res.adj(), res.val()), 1 - res.val()),
            3.0 * 0.07056 * square(A.val()) + 1.5976);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

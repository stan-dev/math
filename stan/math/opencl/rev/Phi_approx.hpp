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
        // e / (1 + e)^2, not res * (1 - res), which cancels for u >> 0
        auto x_sq = square(A.val());
        auto e = exp(-fabs(elt_multiply(A.val(), 0.07056 * x_sq + 1.5976)));
        A.adj() += elt_multiply(
            elt_multiply(res.adj(), elt_divide(e, square(1.0 + e))),
            3.0 * 0.07056 * x_sq + 1.5976);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

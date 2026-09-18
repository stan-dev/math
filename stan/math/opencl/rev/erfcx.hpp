#ifndef STAN_MATH_OPENCL_REV_ERFCX_HPP
#define STAN_MATH_OPENCL_REV_ERFCX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `erfcx()` of a var_value<matrix_cl<double>>.
 *
 * The derivative `2 * x * erfcx(x) - 2 / sqrt(pi)` reuses the function
 * value, so no second `exp` or `erfc` evaluation is needed. That difference
 * cancels for `x >= 4`, where both terms approach `2 / sqrt(pi)` while the
 * result decays like `1 / (sqrt(pi) * x^2)`; the device function
 * `erfcx_derivative` takes the derivative from the tail rational there
 * instead. Without that the error reaches 2.55e+11 ulp at `x = 1e6`. The
 * CPU implementation branches at the same point.
 *
 * @param A argument
 * @return Elementwise `erfcx()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> erfcx(const var_value<T>& A) {
  return make_callback_var(
      erfcx(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(res.adj(),
                                erfcx_derivative(A.val(), res.val()));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

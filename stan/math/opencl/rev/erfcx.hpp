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
 * value, so no second `exp` or `erfc` evaluation is needed and the
 * derivative inherits the value's accuracy in the tails.
 *
 * @param A argument
 * @return Elementwise `erfcx()` of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> erfcx(const var_value<T>& A) {
  return make_callback_var(
      erfcx(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(
            res.adj(), elt_multiply(2.0, elt_multiply(A.val(), res.val()))
                           - TWO_OVER_SQRT_PI);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

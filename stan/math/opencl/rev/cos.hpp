#ifndef STAN_MATH_OPENCL_REV_COS_HPP
#define STAN_MATH_OPENCL_REV_COS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `cos()` of a var_value<matrix_cl<double>>
 * in radians.
 *
 * @param A argument
 * @return Elementwise `cos()` of the input, in radians.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> cos(const var_value<T>& A) {
  return make_callback_var(cos(A.val()),
                           [A](vari_value<matrix_cl<double>>& res) mutable {
                             A.adj() -= elt_multiply(res.adj(), sin(A.val()));
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif

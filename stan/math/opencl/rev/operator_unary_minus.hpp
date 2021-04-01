#ifndef STAN_MATH_OPENCL_REV_OPERATOR_UNARY_MINUS_HPP
#define STAN_MATH_OPENCL_REV_OPERATOR_UNARY_MINUS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the unary minus of the input.
 *
 * @param M input kernel expression
 * @return result of unary minus of the input.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> operator-(const var_value<T>& M) {
  return make_callback_var(-M.val(),
                           [M](vari_value<matrix_cl<double>>& res) mutable {
                             M.adj() -= res.adj();
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif

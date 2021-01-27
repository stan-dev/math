#ifndef STAN_MATH_OPENCL_REV_CROSSPROD_HPP
#define STAN_MATH_OPENCL_REV_CROSSPROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/crossprod.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of pre-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> crossprod(const var_value<T>& M) {
  return make_callback_var(transpose(M.val()) * M.val(),
                           [M](vari_value<matrix_cl<double>>& res) mutable {
                             M.adj() += M.val()
                                        * (res.adj() + transpose(res.adj()));
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif

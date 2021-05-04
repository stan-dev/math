#ifndef STAN_MATH_OPENCL_REV_DIAG_MATRIX_HPP
#define STAN_MATH_OPENCL_REV_DIAG_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return a square diagonal matrix with the specified vector of
 * coefficients as the diagonal values.
 *
 * @param[in] v Specified vector.
 * @return Diagonal matrix with vector as diagonal values.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> diag_matrix(const var_value<T>& v) {
  return make_callback_var(diag_matrix(v.val()),
                           [v](vari_value<matrix_cl<double>>& res) mutable {
                             v.adj() += diagonal(res.adj());
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif

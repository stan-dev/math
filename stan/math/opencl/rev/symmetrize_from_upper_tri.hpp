#ifndef STAN_MATH_OPENCL_REV_SYMMETRIZE_FROM_UPPER_TRI_HPP
#define STAN_MATH_OPENCL_REV_SYMMETRIZE_FROM_UPPER_TRI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return a symmetric matrix using elements from the upper triangular part of
 * the input matrix.
 *
 * @tparam T_x type of elements in the matrix
 * @param A Matrix.
 * @throw std:invalid_argument if the matrix is not square.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> symmetrize_from_upper_tri(
    const var_value<T>& A) {
  return make_callback_var(
      symmetrize_from_upper_tri(A.val()),
      [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += select(row_index() < col_index(),
                          res.adj() + transpose(res.adj()),
                          select(row_index() == col_index(), res.adj(), 0.0));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif

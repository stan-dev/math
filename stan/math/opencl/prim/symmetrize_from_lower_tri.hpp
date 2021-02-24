#ifndef STAN_MATH_OPENCL_PRIM_SYMMETRIZE_FROM_LOWER_TRI_HPP
#define STAN_MATH_OPENCL_PRIM_SYMMETRIZE_FROM_LOWER_TRI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return a symmetric matrix using elements from the upper triangular part of
 * the input matrix.
 *
 * @tparam T_x type of elements in the matrix
 * @param x Matrix.
 * @throw std:invalid_argument if the matrix is not square.
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto symmetrize_from_lower_tri(T_x&& x) {
  check_square("symmetrize_from_lower_tri", "x", x);
  return make_holder_cl(
      [](auto& arg) {
        return select(row_index() < col_index(), transpose(arg), arg);
      },
      std::forward<T_x>(x));
}
}  // namespace math
}  // namespace stan
#endif
#endif

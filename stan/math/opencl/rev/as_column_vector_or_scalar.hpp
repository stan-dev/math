#ifndef STAN_MATH_OPENCL_REV_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_OPENCL_REV_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return a nrows x ncols submatrix starting at (i-1, j-1).
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @param i Starting row.
 * @param j Starting column.
 * @param nrows Number of rows in as_column_vector_or_scalar.
 * @param ncols Number of columns in as_column_vector_or_scalar.
 * @throw std::out_of_range if either index is out of range.
 */
template <typename T,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_any_var_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(const T& m) {
  return m.as_column_vector_or_scalar();
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_PRIM_MAT_FUN_CONCATENATE_ROW_HPP
#define STAN_MATH_PRIM_MAT_FUN_CONCATENATE_ROW_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dims.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Return the result of stacking all elements in the array
 * row-wise. All elements in the array must have matching number of
 * columns.
 *
 * @tparam T Scalar type.
 * @tparam R Row specification.
 * @tparam C Column specification.
 * @param M Array of all elements to be stacked.
 * @return Result of stacking all elements of the array row-wise.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, Eigen::Dynamic, C> concatenate_row(
    const std::vector<Eigen::Matrix<T, R, C>>& M) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  const std::vector<int> M_dims = dims(M);
  std::vector<int> M_rows(M_dims[0]);

  for (int i = 0; i < M_dims[0]; ++i) {
    M_rows[i] = M[i].rows();
    check_size_match("append_row", "columns of M[0]", M_dims[2],
                     "columns of all elements", M[i].cols());
  }

  const int num_rows = std::accumulate(M_rows.begin(), M_rows.end(), 0);

  Matrix<T, Dynamic, C> result(num_rows, M_dims[2]);
  for (int i = 0, offset = 0; i < M_dims[0]; offset += M_rows[i], ++i) {
    result.block(offset, 0, M_rows[i], M_dims[2]) = M[i];
  }
  return result;
}

/**
 * Return the result of stacking all elements in the array
 * row-wise. Whenever all elements are a scalar then this is the same
 * as converting to a vector.
 *
 * @tparam T Scalar type.
 * @param M Array of all elements to be stacked.
 * @return Result of stacking all elements of the array row-wise.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> concatenate_row(
    const std::vector<T>& M) {
  return to_vector(M);
}

}  // namespace math
}  // namespace stan

#endif

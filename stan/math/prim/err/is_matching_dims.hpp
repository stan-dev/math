#ifndef STAN_MATH_PRIM_ERR_IS_MATCHING_DIMS_HPP
#define STAN_MATH_PRIM_ERR_IS_MATCHING_DIMS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/is_size_match.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the two matrices are of the same size.
 * This function checks the runtime sizes only.
 * @tparam T1 Scalar type of the first matrix, requires class method
 *   <code>.size()</code>
 * @tparam T2 Scalar type of the second matrix, requires class method
 *   <code>.size()</code>
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second matrix, can be Eigen::Dynamic
 * @param y1 first matrix to test
 * @param y2 second matrix to test
 * @return <code>true</code> if the dimensions of the matrices match
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline bool is_matching_dims(const Eigen::Matrix<T1, R1, C1>& y1,
                             const Eigen::Matrix<T2, R2, C2>& y2) {
  return is_size_match(y1.rows(), y2.rows())
         && is_size_match(y1.cols(), y2.cols());
}

/**
 * Return <code>true</code> if the two matrices are of the same size.
 * This function checks the runtime sizes and can also check the static
 * sizes as well. For example, a 4x1 matrix is not the same as a vector
 * with 4 elements.
 * @tparam check_compile Whether to check the static sizes
 * @tparam T1 Scalar type of the first matrix
 * @tparam T2 Scalar type of the second matrix
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second matrix, can be Eigen::Dynamic
 * @param y1 first matrix to test
 * @param y2 second matrix to test
 * @return <code>true</code> if the dimensions of the matrices match
 */
template <bool check_compile, typename T1, typename T2, int R1, int C1, int R2,
          int C2>
inline bool is_matching_dims(const Eigen::Matrix<T1, R1, C1>& y1,
                             const Eigen::Matrix<T2, R2, C2>& y2) {
  return !(check_compile && (R1 != R2 || C1 != C2)) && is_matching_dims(y1, y2);
}

}  // namespace math
}  // namespace stan
#endif

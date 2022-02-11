#ifndef STAN_MATH_PRIM_FUN_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_FUN_ROWS_DOT_PRODUCT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of rows of the specified matrices.
 *
 * @tparam Mat1 type of the first matrix (must be derived from \c
 * Eigen::MatrixBase)
 * @tparam Mat2 type of the second matrix (must be derived from \c
 * Eigen::MatrixBase)
 *
 * @param v1 Matrix of first vectors.
 * @param v2 Matrix of second vectors.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the matrices are not the same
 * size
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline Eigen::Matrix<return_type_t<Mat1, Mat2>, Mat1::RowsAtCompileTime, 1>
rows_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("rows_dot_product", "v1", v1, "v2", v2);
  return (v1.cwiseProduct(v2)).rowwise().sum();
}

}  // namespace math
}  // namespace stan

#endif

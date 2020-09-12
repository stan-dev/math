#ifndef STAN_MATH_PRIM_FUN_TO_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_TO_ROW_VECTOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

// row_vector to_row_vector(matrix)
// row_vector to_row_vector(vector)
// row_vector to_row_vector(row_vector)
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, 1, Eigen::Dynamic> to_row_vector(
    const EigMat& matrix) {
  using T = value_type_t<EigMat>;
  Eigen::Matrix<T, 1, Eigen::Dynamic> res(matrix.size());
  Eigen::Map<
      Eigen::Matrix<T, EigMat::RowsAtCompileTime, EigMat::ColsAtCompileTime>>
      res_map(res.data(), matrix.rows(), matrix.cols());
  res_map = matrix;
  return res;
}

// row_vector to_row_vector(real[])
template <typename T>
inline Eigen::Matrix<T, 1, Eigen::Dynamic> to_row_vector(
    const std::vector<T>& vec) {
  return Eigen::Matrix<T, 1, Eigen::Dynamic>::Map(vec.data(), vec.size());
}

// row_vector to_row_vector(int[])
inline Eigen::Matrix<double, 1, Eigen::Dynamic> to_row_vector(
    const std::vector<int>& vec) {
  int C = vec.size();
  Eigen::Matrix<double, 1, Eigen::Dynamic> result(C);
  for (int i = 0; i < C; i++) {
    result(i) = vec[i];
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif

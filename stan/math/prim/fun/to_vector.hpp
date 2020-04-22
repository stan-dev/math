#ifndef STAN_MATH_PRIM_FUN_TO_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_TO_VECTOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

// vector to_vector(matrix)
// vector to_vector(row_vector)
// vector to_vector(vector)
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, 1> to_vector(
    const EigMat& matrix) {
  using T = value_type_t<EigMat>;
  Eigen::Matrix<T, Eigen::Dynamic, 1> res(matrix.size());
  Eigen::Map<
      Eigen::Matrix<T, EigMat::RowsAtCompileTime, EigMat::ColsAtCompileTime>>
      res_map(res.data(), matrix.rows(), matrix.cols());
  res_map = matrix;
  return res;
}

// vector to_vector(real[])
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> to_vector(
    const std::vector<T>& vec) {
  return Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(vec.data(), vec.size());
}

// vector to_vector(int[])
inline Eigen::Matrix<double, Eigen::Dynamic, 1> to_vector(
    const std::vector<int>& vec) {
  int R = vec.size();
  Eigen::Matrix<double, Eigen::Dynamic, 1> result(R);
  for (int i = 0; i < R; i++) {
    result(i) = vec[i];
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif

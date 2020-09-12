#ifndef STAN_MATH_PRIM_FUN_TO_ARRAY_1D_HPP
#define STAN_MATH_PRIM_FUN_TO_ARRAY_1D_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

// real[] to_array_1d(matrix)
// real[] to_array_1d(row_vector)
// real[] to_array_1d(vector)
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline std::vector<value_type_t<EigMat>> to_array_1d(const EigMat& matrix) {
  using T_val = value_type_t<EigMat>;
  std::vector<T_val> result(matrix.size());
  Eigen::Map<Eigen::Matrix<T_val, EigMat::RowsAtCompileTime,
                           EigMat::ColsAtCompileTime>>(
      result.data(), matrix.rows(), matrix.cols())
      = matrix;
  return result;
}

// real[] to_array_1d(...)
template <typename T>
inline std::vector<T> to_array_1d(const std::vector<T>& x) {
  return x;
}

// real[] to_array_1d(...)
template <typename T>
inline std::vector<typename scalar_type<T>::type> to_array_1d(
    const std::vector<std::vector<T>>& x) {
  size_t size1 = x.size();
  size_t size2 = 0;
  if (size1 != 0) {
    size2 = x[0].size();
  }
  std::vector<T> y(size1 * size2);
  for (size_t i = 0, ij = 0; i < size1; i++) {
    for (size_t j = 0; j < size2; j++, ij++) {
      y[ij] = x[i][j];
    }
  }
  return to_array_1d(y);
}

}  // namespace math
}  // namespace stan
#endif

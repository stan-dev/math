#ifndef STAN_MATH_PRIM_FUN_TO_ARRAY_2D_HPP
#define STAN_MATH_PRIM_FUN_TO_ARRAY_2D_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

// real[, ] to_array_2d(matrix)
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline std::vector<std::vector<value_type_t<EigMat>>> to_array_2d(
    const EigMat& matrix) {
  using std::vector;
  using T = value_type_t<EigMat>;
  const Eigen::Ref<const Eigen::Matrix<T, EigMat::RowsAtCompileTime,
                                       EigMat::ColsAtCompileTime>>& mat_ref
      = matrix;
  int C = matrix.cols();
  int R = matrix.rows();
  vector<vector<T>> result(R, vector<T>(C));
  for (int i = 0, ij = 0; i < C; i++) {
    for (int j = 0; j < R; j++, ij++) {
      result[j][i] = mat_ref.data()[ij];
    }
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_MAT_FUN_CONCATENATE_ROW_HPP
#define STAN_MATH_PRIM_MAT_FUN_CONCATENATE_ROW_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dims.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

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

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> concatenate_row(
    const std::vector<T>& M) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  Matrix<T, Dynamic, 1> result(M.size(), 1);
  for (int i = 0; i < M.size(); ++i) {
    result(i) = M[i];
  }
  return result;
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_FUN_ADD_DIAG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element(s) to add along the diagonal
 *
 * @param mat a matrix
 * @param to_add scalar value or column vector or row vector to add along the
 * diagonal
 * @return a matrix with to_add added along main diagonal
 */
template <typename T_m, typename T_a, typename = require_eigen_t<T_m>,
          typename = require_any_t<is_eigen_vector<T_a>, is_stan_scalar<T_a>>>
inline typename Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic,
                              Eigen::Dynamic>
add_diag(const T_m &mat, const T_a &to_add) {
  if (is_vector<T_a>::value) {
    const size_t length_diag = std::min(mat.rows(), mat.cols());
    check_consistent_size("add_diag", "number of elements of to_add", to_add,
                          length_diag);
  }
  Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic, Eigen::Dynamic> out
      = mat;
  out.diagonal().array()
      += as_array_or_scalar(as_column_vector_or_scalar(to_add));
  return out;
}

}  // namespace math
}  // namespace stan

#endif

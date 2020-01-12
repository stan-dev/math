#ifndef STAN_MATH_PRIM_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_FUN_ADD_DIAG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element to add along the diagonal
 *
 * @param mat a matrix
 * @param to_add value to add along the diagonal
 * @return a matrix with to_add added along main diagonal
 */
template <typename T_m, typename T_a>
inline typename Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic,
                              Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, Eigen::Dynamic, Eigen::Dynamic> &mat,
         const T_a &to_add) {
  Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic, Eigen::Dynamic> out
      = mat;
  out.diagonal().array() += to_add;
  return out;
}

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element to add along the diagonal
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param mat a matrix
 * @param to_add Sequence of values to add along the diagonal
 * @return a matrix with to_add added along main diagonal
 * @throw invalid_argument if to_add is vector-like but does not have
 * the same number of elements as the main diagonal of mat
 */
template <typename T_m, typename T_a, int R, int C>
inline typename Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic,
                              Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, Eigen::Dynamic, Eigen::Dynamic> &mat,
         const Eigen::Matrix<T_a, R, C> &to_add) {
  const size_t length_diag = std::min(mat.rows(), mat.cols());
  check_consistent_size("add_diag", "number of elements of to_add", to_add,
                        length_diag);

  Eigen::Matrix<return_type_t<T_m, T_a>, Eigen::Dynamic, Eigen::Dynamic> out
      = mat;
  out.diagonal() += to_add;
  return out;
}

}  // namespace math
}  // namespace stan

#endif

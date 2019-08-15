#ifndef STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <algorithm>

namespace stan {
namespace math {
/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element to add along the diagonal
 * @param mat a matrix
 * @param to_add value to add along the diagonal
 * @return a matrix with to_add added along main diagonal
 */
template <typename T_m, typename T_a, enable_if_eigen<T_m> * = nullptr,
          enable_if_not_eigen<T_a> * = nullptr>
inline auto add_diag(const T_m &mat, const T_a &to_add) {
  T_m out = mat;
  out.diagonal().array() += to_add;
  return out;
}

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element to add along the diagonal
 * @param mat a matrix
 * @param to_add Sequence of values to add along the diagonal
 * @return a matrix with to_add added along main diagonal
 * @throw invalid_argument if to_add is vector-like but does not have
 * the same number of elements as the main diagonal of mat
 */
template <typename T_m, typename T_a, enable_if_all_eigen<T_m, T_a> * = nullptr>
inline auto add_diag(const T_m &mat, const T_a &to_add) {
  const size_t length_diag = std::min(mat.rows(), mat.cols());
  check_consistent_size("add_diag", "number of elements of to_add", to_add,
                        length_diag);

  T_m out = mat;
  out.diagonal() += to_add;
  return out;
}
}  // namespace math
}  // namespace stan
#endif

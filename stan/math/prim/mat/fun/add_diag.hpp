#ifndef STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <algorithm>

namespace stan {
namespace math {
/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a Type of element to add along the diagonal
 * @param mat a matrix
 * @param to_add (Sequence of) value(s) to add along the diagonal
 * @return a matrix with to_add added along main diagonal
 * @throw invalid_argument if to_add is vector-like but does not have
 * the same number of elements as the main diagonal of mat
 */
template <typename T_m, typename T_a>
inline typename Eigen::Matrix<typename return_type<T_m, T_a>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, Eigen::Dynamic, Eigen::Dynamic> &mat,
         const T_a &to_add) {
  size_t length_diag = std::min(mat.rows(), mat.cols());
  check_consistent_size("add_diag", "number of elements of to_add", to_add,
                        length_diag);

  Eigen::Matrix<typename return_type<T_m, T_a>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      out = mat;
  scalar_seq_view<T_a> to_add_vec(to_add);

  for (size_t i = 0; i < length_diag; ++i)
    out(i, i) += to_add_vec[i];

  return out;
}
}  // namespace math
}  // namespace stan
#endif

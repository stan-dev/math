#ifndef STAN_MATH_PRIM_FUN_TO_ROW_VECTOR_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_TO_ROW_VECTOR_ARRAY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a standard vector of Eigen row vectors from the rows of the input
 * matrix.
 *
 * @tparam EigMat type of the input matrix
 * @param matrix input matrix
 * @return the array of row vectors representation of the input
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline std::vector<Eigen::Matrix<value_type_t<EigMat>, 1, Eigen::Dynamic>>
to_row_vector_array(const EigMat& matrix) {
  using T = value_type_t<EigMat>;
  std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic>> result;
  result.reserve(matrix.rows());
  for (int i = 0; i < matrix.rows(); ++i) {
    result.push_back(matrix.row(i));
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif

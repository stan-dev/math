#ifndef STAN_MATH_PRIM_FUN_TO_VECTOR_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_TO_VECTOR_ARRAY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a standard vector of Eigen column vectors from the columns of
 * the input matrix.
 *
 * @tparam EigMat type of the input matrix
 * @param matrix input matrix
 * @return the array of column vectors representation of the input
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline std::vector<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, 1>>
to_vector_array(const EigMat& matrix) {
  using T = value_type_t<EigMat>;
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> result;
  result.reserve(matrix.cols());
  for (int j = 0; j < matrix.cols(); ++j) {
    result.push_back(matrix.col(j));
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif

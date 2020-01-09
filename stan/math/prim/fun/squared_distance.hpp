#ifndef STAN_MATH_PRIM_MAT_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <int R, int C>
inline double squared_distance(const Eigen::Matrix<double, R, C>& v1,
                               const Eigen::Matrix<double, R, C>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return (v1 - v2).squaredNorm();
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam R1 number of rows in the first vector, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first vector, can be Eigen::Dynamic
 * @tparam R2 number of rows in the second vector, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second vector, can be Eigen::Dynamic
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <int R1, int C1, int R2, int C2>
inline double squared_distance(const Eigen::Matrix<double, R1, C1>& v1,
                               const Eigen::Matrix<double, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return (v1.transpose() - v2).squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif

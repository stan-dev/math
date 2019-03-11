#ifndef STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 *
 * @param x1 First vector.
 * @param x2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T1, typename T2>
inline typename boost::math::tools::promote_args<T1, T2>::type squared_distance(
    const T1& x1, const T2& x2) {
  check_finite("squared_distance", "x1", x1);
  check_finite("squared_distance", "x2", x2);
  return square(x1 - x2);
}
}  // namespace math
}  // namespace stan

namespace stan {
namespace math {

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam R Rows at compile time of vector inputs
 * @tparam C columns at compile time of vector inputs
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
 * @tparam R1 Rows at compile time of first vector input
 * @tparam C1 Columns at compile time of first vector input
 * @tparam R2 Rows at compile time of second vector input
 * @tparam C2 Columns at compile time of second vector input
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

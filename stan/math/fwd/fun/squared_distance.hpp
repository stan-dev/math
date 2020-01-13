#ifndef STAN_MATH_FWD_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_FWD_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/fwd/fun/dot_self.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T, int R, int C>
inline fvar<T> squared_distance(const Eigen::Matrix<fvar<T>, R, C>& v1,
                                const Eigen::Matrix<double, R, C>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, R, C> v3 = subtract(v1, v2);
  return dot_self(v3);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
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
template <typename T, int R1, int C1, int R2, int C2>
inline fvar<T> squared_distance(const Eigen::Matrix<fvar<T>, R1, C1>& v1,
                                const Eigen::Matrix<double, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<double, R1, C1> t_v2 = v2.transpose();
  Eigen::Matrix<fvar<T>, R1, C1> v3 = subtract(v1, t_v2);
  return dot_self(v3);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T, int R, int C>
inline fvar<T> squared_distance(const Eigen::Matrix<double, R, C>& v1,
                                const Eigen::Matrix<fvar<T>, R, C>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, R, C> v3 = subtract(v1, v2);
  return dot_self(v3);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
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
template <typename T, int R1, int C1, int R2, int C2>
inline fvar<T> squared_distance(const Eigen::Matrix<double, R1, C1>& v1,
                                const Eigen::Matrix<fvar<T>, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<double, R2, C2> t_v1 = v1.transpose();
  Eigen::Matrix<fvar<T>, R2, C2> v3 = subtract(t_v1, v2);
  return dot_self(v3);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T, int R, int C>
inline fvar<T> squared_distance(const Eigen::Matrix<fvar<T>, R, C>& v1,
                                const Eigen::Matrix<fvar<T>, R, C>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, R, C> v3 = subtract(v1, v2);
  return dot_self(v3);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T inner type of the fvar vector
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
template <typename T, int R1, int C1, int R2, int C2>
inline fvar<T> squared_distance(const Eigen::Matrix<fvar<T>, R1, C1>& v1,
                                const Eigen::Matrix<fvar<T>, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, R2, C2> t_v1 = v1.transpose();
  Eigen::Matrix<fvar<T>, R2, C2> v3 = subtract(t_v1, v2);
  return dot_self(v3);
}

}  // namespace math
}  // namespace stan
#endif

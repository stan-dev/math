#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam T type of elements in the vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 */
template <typename T, int R, int C>
Eigen::Matrix<T, R, C> unit_vector_constrain(const Eigen::Matrix<T, R, C>& y) {
  using std::sqrt;
  check_vector("unit_vector_constrain", "y", y);
  check_nonzero_size("unit_vector_constrain", "y", y);
  T SN = dot_self(y);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  return y / sqrt(SN);
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam T type of elements in the vector
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @param lp Log probability reference to increment.
 */
template <typename T, int R, int C>
Eigen::Matrix<T, R, C> unit_vector_constrain(const Eigen::Matrix<T, R, C>& y,
                                             T& lp) {
  using std::sqrt;
  check_vector("unit_vector_constrain", "y", y);
  check_nonzero_size("unit_vector_constrain", "y", y);
  T SN = dot_self(y);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  lp -= 0.5 * SN;
  return y / sqrt(SN);
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_MEAN_HPP
#define STAN_MATH_PRIM_FUN_MEAN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample mean (i.e., average) of the coefficients
 * in the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param v Specified vector.
 * @return Sample mean of vector coefficients.
 * @throws std::domain_error if the size of the vector is less
 * than 1.
 */
template <typename T>
inline return_type_t<T> mean(const std::vector<T>& v) {
  check_nonzero_size("mean", "v", v);
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> m(&v[0], v.size());
  return m.mean();
}

/**
 * Returns the sample mean (i.e., average) of the coefficients
 * in the specified vector, row vector, or matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Specified vector, row vector, or matrix.
 * @return Sample mean of vector coefficients.
 */
template <typename T, int R, int C>
inline return_type_t<T> mean(const Eigen::Matrix<T, R, C>& m) {
  check_nonzero_size("mean", "m", m);
  return m.mean();
}

}  // namespace math
}  // namespace stan

#endif

#include <vector>

#ifndef STAN_MATH_PRIM_MAT_META_DUPLICATE_IF_SCALAR_HPP
#define STAN_MATH_PRIM_MAT_META_DUPLICATE_IF_SCALAR_HPP

namespace stan {
namespace math {
/**
 * This program is used to either duplicate a scalar to a vector of a given
 * length or to simply act as the identity on a vector input.
 *
 * @tparam T Type of scalar entries.
 * @param arg1 Argument to be duplicated.
 * @param N Number of duplications to be made.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
  duplicate_if_scalar(const T &arg1, int N) {
    return Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(N, 1) * arg1;
  }

/**
 * This program is used to either duplicate a scalar to a vector of a given
 * length or to simply act as the identity on a vector input.
 *
 * @tparam T Type of vector entries.
 * @param arg1 Argument to be returned.
 * @param N Number of duplications to be made (ignored).
 */
template <typename T>
const Eigen::Matrix<T, Eigen::Dynamic, 1>&
  duplicate_if_scalar(const Eigen::Matrix<T, Eigen::Dynamic, 1> &arg1, int N) {
    return arg1;
  }
}  // namespace math
}  // namespace stan

#endif

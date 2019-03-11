#ifndef STAN_MATH_PRIM_FUN_SUM_HPP
#define STAN_MATH_PRIM_FUN_SUM_HPP

#include <cstddef>
#include <vector>
#include <numeric>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sum.hpp>

namespace stan {
namespace math {

/**
 * Return the sum of the values in the specified standard vector.
 *
 * @tparam T Type of elements summed.
 * @param xs Standard vector to sum.
 * @return Sum of elements.
 */
template <typename T>
inline T sum(const std::vector<T>& xs) {
  return std::accumulate(xs.begin(), xs.end(), T{0});
}

}  // namespace math
}  // namespace stan

namespace stan {
namespace math {

/**
 * Returns the sum of the coefficients of the specified
 * column vector.
 *
 * @tparam T Type of elements in matrix.
 * @tparam R Row type of matrix.
 * @tparam C Column type of matrix.
 * @param v Specified vector.
 * @return Sum of coefficients of vector.
 */
template <typename T, int R, int C>
inline T sum(const Eigen::Matrix<T, R, C>& v) {
  return v.sum();
}

}  // namespace math
}  // namespace stan
#endif

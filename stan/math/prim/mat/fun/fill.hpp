#ifndef STAN_MATH_PRIM_MAT_FUN_FILL_HPP
#define STAN_MATH_PRIM_MAT_FUN_FILL_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * The specified matrix is filled by element.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @tparam S Type of value.
 *
 * @param x Container.
 * @param y Value.
 */
template <typename T, int R, int C, typename S>
void fill(Eigen::Matrix<T, R, C>& x, const S& y) {
  x.fill(y);
}

}  // namespace math
}  // namespace stan

#endif

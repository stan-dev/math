#ifndef STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP
#define STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns 1, the number of elements in a primitive type.
 *
 * @tparam T scalar type
 * @param x Argument of primitive type.
 * @return 1
 */
template <typename T>
inline int num_elements(const T& x) {
  return 1;
}

/**
 * Returns the size of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m argument matrix
 * @return size of matrix
 */
template <typename T, int R, int C>
inline int num_elements(const Eigen::Matrix<T, R, C>& m) {
  return m.size();
}

/**
 * Returns the number of elements in the specified vector.
 * This assumes it is not ragged and that each of its contained
 * elements has the same number of elements.
 *
 * @tparam T type of elements in the vector
 * @param v argument vector
 * @return number of contained arguments
 */
template <typename T>
inline int num_elements(const std::vector<T>& v) {
  if (v.size() == 0) {
    return 0;
  }
  return v.size() * num_elements(v[0]);
}

}  // namespace math
}  // namespace stan

#endif

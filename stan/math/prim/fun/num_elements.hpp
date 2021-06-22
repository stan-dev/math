#ifndef STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP
#define STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
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
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline int num_elements(const T& x) {
  return 1;
}

/**
 * Returns the size of the specified matrix.
 *
 * @tparam T type of the matrix
 *
 * @param m argument matrix
 * @return size of matrix
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline int num_elements(const T& m) {
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

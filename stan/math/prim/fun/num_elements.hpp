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
 * @return 1
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline constexpr int num_elements(const T& /* x*/) {
  return 1;
}

/**
 * Returns the size of the specified matrix.
 *
 * @tparam T type of the matrix
 * @param x argument matrix
 * @return size of matrix
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline int num_elements(const T& x) {
  return x.size();
}

/**
 * Returns the number of elements in the specified vector.
 * This assumes it is not ragged and that each of its contained
 * elements has the same number of elements.
 *
 * @tparam T type of elements in the vector
 * @param x argument vector
 * @return number of contained arguments
 */
template <typename T>
inline int num_elements(const std::vector<T>& x) {
  if (x.size() == 0) {
    return 0;
  }
  return x.size() * num_elements(x[0]);
}

}  // namespace math
}  // namespace stan

#endif

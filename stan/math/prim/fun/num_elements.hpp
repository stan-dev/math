#ifndef STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP
#define STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
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
inline size_t num_elements(const T& x) {
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
inline size_t num_elements(const T& m) {
  return m.size();
}

/**
 * Returns the number of elements in the specified vector.
 * @tparam T type of elements in the vector
 * @param v argument vector
 * @return number of contained arguments
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline size_t num_elements(const std::vector<T>& v) {
  return v.size();
}

/**
 * Returns the number of elements in the specified vector
 *
 * @tparam T type of elements in the vector
 * @param v argument vector
 * @return number of contained arguments
 */
template <typename T, require_container_t<T>* = nullptr>
inline size_t num_elements(const std::vector<T>& v) {
  size_t size = 0;
  for (auto&& v_val : v) {
    size += num_elements(v_val);
  }
  return size;
}

/**
 * Returns the number of elements in the specified tuple
 *
 * @tparam T type of tuple
 * @param v tuple
 * @return number of contained arguments
 */
template <typename T, require_tuple_t<T>* = nullptr>
inline size_t num_elements(const T& v) {
  size_t size = 0;
  math::apply(
      [&size](const auto&... args) {
        static_cast<void>(
            std::initializer_list<int>{(size += num_elements(args), 0)...});
      },
      v);
  return size;
}

}  // namespace math
}  // namespace stan

#endif

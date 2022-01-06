#ifndef STAN_MATH_PRIM_FUN_MAX_HPP
#define STAN_MATH_PRIM_FUN_MAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the maximum value of the two specified
 * scalar arguments.
 *
 * @tparam T1 type of first argument (must be arithmetic)
 * @tparam T2 type of second argument (must be arithmetic)
 * @param x first argument
 * @param y second argument
 * @return maximum value of the two arguments
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
return_type_t<T1, T2> max(T1 x, T2 y) {
  return std::max(x, y);
}

/**
 * Returns the maximum value in the specified
 * matrix, vector, row vector or std vector.
 *
 * @tparam T type of the container
 * @param m specified matrix, vector, row vector or std vector
 * @return maximum value in the container, or -infinity if the
 * container is size zero and the scalar type in container is floating point
 * number
 * @throws <code>std::invalid_argument</code> if the vector is size zero and the
 * scalar type in the container is integer
 */
template <typename T, require_container_t<T>* = nullptr>
inline value_type_t<T> max(const T& m) {
  if (std::is_integral<value_type_t<T>>::value) {
    check_nonzero_size("max", "int vector", m);
  } else if (m.size() == 0) {
    return NEGATIVE_INFTY;
  }
  return apply_vector_unary<T>::reduce(
      m, [](const auto& x) { return x.maxCoeff(); });
}

}  // namespace math
}  // namespace stan

#endif

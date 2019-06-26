#ifndef STAN_MATH_PRIM_ARR_FUN_VEC_CONCAT_HPP
#define STAN_MATH_PRIM_ARR_FUN_VEC_CONCAT_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Ends the recursion to extract the event stack.
 * @param v1 Events on the OpenCL event stack.
 * @return returns input vector.
 */
template <typename T>
inline const std::vector<T>& vec_concat(const std::vector<T>& v1) {
  return v1;
}

/**
 * Gets the event stack from a vector of events and other arguments.
 * @param v1 A event stack to roll up.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam Args Types for variadic.
 * @return Vector of OpenCL events
 */
template <typename T, typename... Args>
inline const std::vector<T> vec_concat(const std::vector<T>& v1,
                                       const Args... args) {
  std::vector<T> vec = vec_concat(args...);
  vec.insert(vec.end(), v1.begin(), v1.end());
  return vec;
}
}  // namespace math
}  // namespace stan

#endif

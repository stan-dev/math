#ifndef STAN_MATH_PRIM_FUN_VEC_CONCAT_HPP
#define STAN_MATH_PRIM_FUN_VEC_CONCAT_HPP

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
template <typename Vec>
inline decltype(auto) vec_concat(Vec&& v1) {
  return std::forward<Vec>(v1);
}

/**
 * Get the event stack from a vector of events and other arguments.
 *
 * @tparam Vec type of first input standard vector
 * @tparam VecArgs types for standard vector variadic arguments
 * @param v1 event stack to roll up
 * @param args variadic arguments passed down to the next recursion
 * @return Vector of OpenCL events
 */
template <typename Vec, typename... VecArgs, require_all_std_vector_t<Vec, VecArgs...>* = nullptr>
inline auto vec_concat(Vec&& v1, VecArgs&&... args) {
  std::vector<value_type_t<Vec>> vec = vec_concat(std::forward<VecArgs>(args)...);
  vec.insert(vec.end(), v1.begin(), v1.end());
  return vec;
}

}  // namespace math
}  // namespace stan

#endif

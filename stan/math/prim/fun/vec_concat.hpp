#ifndef STAN_MATH_PRIM_FUN_VEC_CONCAT_HPP
#define STAN_MATH_PRIM_FUN_VEC_CONCAT_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {
namespace internal {
inline auto sum_vector_sizes() { return 0; }
/**
 * Get the internal sizes of a pack of vectors.
 * @tparam Vec Type of first vector to count.
 * @tparam VecArgs Parameter pack of vector types to accumulate sizes of
 * @param x vector to start accumulation of sizes with.
 * @param args pack of vectors to accumulate sizes for.
 */
template <typename Vec, typename... VecArgs>
inline auto sum_vector_sizes(const Vec& x, const VecArgs&... args) {
  return x.size() + sum_vector_sizes(args...);
}
/**
 * End of recursion for appending to vector.
 */
template <typename VecInOut>
inline void append_vectors(VecInOut& x) {}
/**
 * Fill a vector with other vectors.
 * @tparam VecInOut Type of vector to be filled.
 * @tparam VecIn Type of the first vector to insert into the in/out vector.
 * @tparam VecArgs Parameter pack of other vectors to fill in the in/out vector.
 * @param x Vector to be filled.
 * @param y First vector to insert into the in/out vector.
 * @param args Pack of other vectors to fill in the in/out vector.
 */
template <typename VecInOut, typename VecIn, typename... VecArgs>
inline void append_vectors(VecInOut& x, const VecIn& y,
                           const VecArgs&... args) {
  x.insert(x.end(), y.begin(), y.end());
  append_vectors(x, args...);
}
}  // namespace internal

/**
 * Get the event stack from a vector of events and other arguments.
 *
 * @tparam Vec type of elements in the array
 * @tparam Args types for variadic arguments
 * @param v1 event stack to roll up
 * @param args variadic arguments passed down to the next recursion
 * @return Vector of OpenCL events
 */
template <typename Vec, typename... Args>
inline auto vec_concat(const Vec& v1, const Args&... args) {
  std::vector<value_type_t<Vec>> vec;
  vec.reserve(internal::sum_vector_sizes(v1, args...));
  internal::append_vectors(vec, v1, args...);
  return vec;
}

}  // namespace math
}  // namespace stan

#endif

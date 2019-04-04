#ifndef STAN_MATH_PRIM_ARR_FUN_VEC_CONCAT_HPP
#define STAN_MATH_PRIM_ARR_FUN_VEC_CONCAT_HPP

#include <type_traits>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename T>
void concat_helper(std::vector<T>& l, const std::vector<T>& r) {
  l.insert(l.end(), r.begin(), r.end());
}

template <typename T>
void concat_helper(std::vector<T>& l, std::vector<T>&& r) {
  l.insert(l.end(), std::make_move_iterator(r.begin()),
           std::make_move_iterator(r.end()));
}
}  // namespace internal

/**
 * Gets the event stack from a vector of events and other arguments.
 * @param v1 A event stack to roll up.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam Args Types for variadic.
 * @return Vector of OpenCL events
 */
template <typename T, typename... Args>
std::vector<T> vec_concat(const std::vector<T> v1, Args&&... args) {
  std::size_t s = v1.size();
  // Hacks to expand the variadic template
  auto dummy1 = {s += args.size()...};
  std::vector<T> vec;
  vec.reserve(s);
  auto dummy2
      = {(internal::concat_helper(vec, std::forward<Args>(args)), 0)...};
  vec.insert(vec.end(), v1.begin(), v1.end());
  return std::move(vec);
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_TAIL_HPP
#define STAN_MATH_PRIM_FUN_TAIL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector or row vector (same as
 * input) from the back of the specified vector or row vector.
 *
 * @tparam T type of the vector
 * @param v Vector input.
 * @param n Size of return.
 * @return The last n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T, require_vector_t<T>* = nullptr>
inline auto tail(T&& v, size_t n) {
  if constexpr (is_std_vector_v<T>) {
    if (n != 0) {
      check_std_vector_index("tail", "n", v, n);
    }
    std::decay_t<T> s(v.end() - n, v.end());
    return s;
  } else {
    if (n != 0) {
      check_vector_index("tail", "n", v, n);
    }
    return make_holder(
        [n](auto&& v_) { return std::forward<decltype(v_)>(v_).tail(n); },
        std::forward<T>(v));
  }
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_SEGMENT_HPP
#define STAN_MATH_PRIM_FUN_SEGMENT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a row/column vector starting
 * from the specified element - 1 of the specified row/column vector.
 *
 * @tparam T type of the vector
 */
template <typename Vec, require_vector_t<Vec>* = nullptr>
inline auto segment(Vec&& v, size_t i, size_t n) {
  check_greater("segment", "n", i, 0.0);
  check_less_or_equal("segment", "n", i, static_cast<size_t>(v.size()));
  if (n != 0) {
    check_greater("segment", "n", i + n - 1, 0.0);
    check_less_or_equal("segment", "n", i + n - 1,
                        static_cast<size_t>(v.size()));
  }
  if constexpr (is_std_vector_v<Vec>) {
    std::decay_t<Vec> s;
    for (size_t j = 0; j < n; ++j) {
      s.push_back(v[i + j - 1]);
    }
    return s;
  } else {
    return make_holder([i, n](auto&& v_) { return v_.segment(i - 1, n); },
                      std::forward<Vec>(v));
  }
}

}  // namespace math
}  // namespace stan

#endif

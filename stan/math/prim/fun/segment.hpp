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
inline auto segment(const Vec& v, size_t i, size_t n) {
  check_greater("segment", "n", i, 0.0);
  check_less_or_equal("segment", "n", i, static_cast<size_t>(v.size()));
  if (n != 0) {
    check_greater("segment", "n", i + n - 1, 0.0);
    check_less_or_equal("segment", "n", i + n - 1,
                        static_cast<size_t>(v.size()));
  }
  return v.segment(i - 1, n);
}

template <typename T>
std::vector<T> segment(const std::vector<T>& sv, size_t i, size_t n) {
  check_greater("segment", "i", i, 0.0);
  check_less_or_equal("segment", "i", i, sv.size());
  if (n != 0) {
    check_greater("segment", "i+n-1", i + n - 1, 0.0);
    check_less_or_equal("segment", "i+n-1", i + n - 1,
                        static_cast<size_t>(sv.size()));
  }
  std::vector<T> s;
  for (size_t j = 0; j < n; ++j) {
    s.push_back(sv[i + j - 1]);
  }
  return s;
}

}  // namespace math
}  // namespace stan

#endif

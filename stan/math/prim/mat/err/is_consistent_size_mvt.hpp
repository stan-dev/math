#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZE_MVT_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZE_MVT_HPP

#include <stan/math/prim/mat/meta/length.hpp>
#include <stan/math/prim/mat/meta/length_mvt.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Check if the dimension of `x` is consistent, which is defined to be
 * <code>expected_size</code> if `x` is a vector of vectors or 1 if `x` is
 * a single vector.
 * @tparam T Type of value
 * @param x Variable to check for consistent size
 * @param expected_size Expected size if `x` is a vector
 * @return `true` if the size is consistent
 */
template <typename T>
inline bool is_consistent_size_mvt(const T& x, size_t expected_size) {
  return length(x) == expected_size
         || !is_vector<
                typename std::remove_reference<decltype(x[0])>::type>::value;
}

}  // namespace math
}  // namespace stan
#endif

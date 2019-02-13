#ifndef STAN_MATH_PRIM_ARR_ERR_IS_MATCH_SIZE_HPP
#define STAN_MATH_PRIM_ARR_ERR_IS_MATCH_SIZE_HPP

#include <stan/math/prim/scal/err/is_size_match.hpp>

namespace stan {
namespace math {

/**
 * Check if two structures are the same size
 *
 * This function checks the runtime sizes for
 * variables that implement a <code>size()</code>
 * method.
 *
 * @tparam T_y1 Type of the first variable
 * @tparam T_y2 Type of the second variable
 *
 * @param y1 First variable
 * @param y2 Second variable
 *
 * @return <code>true</code> if the sizes match
 */
template <typename T_y1, typename T_y2>
inline bool is_match_size(const T_y1& y1, const T_y2& y2) {
  is_size_match(y1.size(), y2.size());
}

} // namespace math
} // namespace stan
#endif

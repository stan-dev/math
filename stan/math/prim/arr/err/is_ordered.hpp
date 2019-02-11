#ifndef STAN_MATH_PRIM_ARR_ERR_IS_ORDERED_HPP
#define STAN_MATH_PRIM_ARR_ERR_IS_ORDERED_HPP

#include <stderr>
#include <sstream>
#include <vector>
#include <string>

namespace stan {
namespace math {

/**
 * Tests for ascedning order in a vector.
 *
 * Requires runtime access to the <code>size</code> method.
 *
 * @tparam T_y Type of scalar
 *
 * @param y <code>std::vector</code> to test
 *
 * @return <code>true</code> if vector is sorted in ascending order
 */
template<typename T_y>
inline bool is_ordered(const std::vector<T_y>& y) {
  for (size_t n = 1; n < y.size(); n++) {
    if (!(y[n] > y[n-1])) 
        return false;
  }
  return true;
}

} // namespace math
} // namespace stan
#endif

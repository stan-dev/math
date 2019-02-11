#ifndef STAN_MATH_PRIM_ARR_FUN_IS_ARR_SIZE_MATCH.HPP
#define STAN_MATH_PRIM_ARR_FUN_IS_ARR_SIZE_MATCH.HPP

namespace stan {
namespace math {

/**
 * Check if two structures are of the same size
 *
 * Requires runtime  access to the <code>size()</code> method.
 *
 * @tparam T_y1 Type of first variable
 * @tparam T_y2 Type of second variable
 *
 * @param y1 struct First structure
 * @param y2 struct Second structure
 *
 * @return <code>true</code> if the inputs are of the same size
 */
template <typename T_y1, typename T_y2>
inline bool is_arr_size_match(const T_y1& y1,
			  const T_y2& y2) {
  if (y1.size() == y2.size()) 
    return true;
  return false;
}

} // namespace math
} // namespace stan
#endif

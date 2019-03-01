#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_NONZERO_SIZE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_NONZERO_SIZE_HPP

namespace stan {
namespace math {

/**
 * Returns <code>true</code> if the container has size zero.
 *
 * @tparam T_y Type of container, requires class method <code>.size()</code>
 *
 * @param y Container to test -- matrix/vector
 *
 * @return <code>true</code> if container has size zero
 */
template <typename T_y>
inline bool is_nonzero_size(const T_y& y) {
  if (y.size() > 0)
    return true;
  return false;
}

}  // namespace math
}  // namespace stan
#endif

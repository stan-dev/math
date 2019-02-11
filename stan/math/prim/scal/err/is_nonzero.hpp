#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_NONZERO_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_NONZERO_HPP

namespace stan {
namespace math {

/**
 * Checks if the structure is of non-zero size.
 *
 * Requirs runtime access to the <code>size()</code> method.
 *
 * @tparam T_y Type of container
 *
 * @param y structure Structure to test. Matrices or vector.
 *
 * @return <code>true</code> if <code>size()</code> returns
 *   greater than 0
 */
template <typename T_y>
inline bool is_nonzero(const T_y& y) {
  if (y.size() > 0)
    return true;
  return false;
}

}  // namespace math
}  // namespace stan
#endif

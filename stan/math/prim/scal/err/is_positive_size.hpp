#ifndef STAN_MATH_PRIM_SCAL_ERR_CHECK_IS_POSITIVE_SIZE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_CHECK_IS_POSITIVE_SIZE_HPP

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>size</code> is positive.
 * @param size Size value to check
 * @return <code>true</code> if <code>size</code> is not zero or negative.
 */
inline bool is_positive_size(int size) { return size > 0; }

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_SCAL_FUN_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_FUN_AS_ARRAY_OR_SCALAR_HPP

namespace stan {
namespace math {

/**
 * Converts a matrix type to an array. For scalar inputs this is identity
 * function.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T>
inline T as_array_or_scalar(const T& v) {
  return v;
}

}  // namespace math
}  // namespace stan

#endif

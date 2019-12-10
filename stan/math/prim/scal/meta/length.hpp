#ifndef STAN_MATH_PRIM_SCAL_META_LENGTH_HPP
#define STAN_MATH_PRIM_SCAL_META_LENGTH_HPP

#include <cstdlib>
#include <vector>

namespace stan {
/** \ingroup type_trait
 * Returns the length of primitive scalar types
 * that are always of length 1.
 */
template <typename T>
size_t length(const T& /*x*/) {
  return 1U;
}

/** \ingroup type_trait
 * Returns the length of the provided std::vector.
 *
 * @param x input vector
 * @tparam T type of the elements in the vector
 * @return the length of the input vector
 */
template <typename T>
size_t length(const std::vector<T>& x) {
  return x.size();
}
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_MAT_ERR_IS_STD_VECTOR_INDEX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_STD_VECTOR_INDEX_HPP

#include <stan/math/prim/scal/err/out_of_range.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified index is a valid
 * <code>std::</code> vector.
 * This check is 1-indexed by default. This behavior can be changed
 * by setting <code>stan::error_index::value</code>.
 * @tparam T Scalar type, requires class method <code>.size()</code>
 * @param y <code>std::vector</code> to test
 * @param i Index to test
 * @return <code>true</code> if the index is not out of range
 */
template <typename T>
inline bool is_std_vector_index(const std::vector<T>& y, int i) {
  return i >= static_cast<int>(stan::error_index::value)
         && i < static_cast<int>(y.size() + stan::error_index::value);
}

}  // namespace math
}  // namespace stan
#endif

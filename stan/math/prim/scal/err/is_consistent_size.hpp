#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZE_HPP

#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <stan/math/prim/scal/meta/size_of.hpp>

namespace stan {
namespace math {

/**
 * Check if the dimension of <code>x</code> is consistent, which is defined to
 * be <code>expected_size</code> if <code>x</code> is a vector or 1 if
 * <code>x</code> is not a vector.
 * @tparam T Type of value
 * @param x Variable to check for consistent size
 * @param expected_size Expected size if <code>x</code> is a vector
 * @return <code>true</code> if the size is consistent
 */
template <typename T>
inline bool is_consistent_size(const T& x, size_t expected_size) {
  return (!is_vector<T>::value
          || (is_vector<T>::value && expected_size == stan::size_of(x)));
}

}  // namespace math
}  // namespace stan
#endif

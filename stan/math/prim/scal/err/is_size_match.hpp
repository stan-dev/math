#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_SIZE_MATCH_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_SIZE_MATCH_HPP

#include <boost/type_traits/common_type.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>

namespace stan {
namespace math {

/**
 * Check if the provided sizes match
 *
 * @tparam T_size1 Type of size 1
 * @tparam T_size2 Type of size 2
 *
 * @param i Size 1
 * @param j Size 2
 *
 * @return <code>true</code> if provided dimensions match
 */
template <typename T_size1, typename T_size2>
inline bool is_size_match(T_size1 i, T_size2 j) {
  if (likely(i == static_cast<T_size1>(j)))
    return true;
  return false;
}

}  // namespace math
}  // namespace stan
#endif


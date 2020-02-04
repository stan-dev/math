#ifndef STAN_MATH_PRIM_META_SIZE_HPP
#define STAN_MATH_PRIM_META_SIZE_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
namespace math {
/** \ingroup type_trait
 * Return 1, which is the length of scalar types.
 *
 * @return 1
 */
template <typename T, typename = require_stan_scalar_t<T>>
inline size_t size(const T& /*x*/) {
  return 1U;
}

/** \ingroup type_trait
 * Return the size of the specified standard vector.
 *
 * @tparam T value type of vector
 * @param[in] m vector whose size is returned
 */
template <typename T>
inline size_t size(const std::vector<T>& m) {
  return m.size();
}

/**
 * Return the size of the specified Eigen matrix, vector, or row
 * vector.
 *
 * @tparam T value type
 * @tparam R row specification
 * @tparam C column specification
 * @param[in] x matrix, vector, or row vector
 * @return size of argument
 */
template <typename T, int R, int C>
inline size_t size(const Eigen::Matrix<T, R, C>& x) {
  return x.size();
}
}  // namespace math
}  // namespace stan
#endif

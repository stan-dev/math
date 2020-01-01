#ifndef STAN_MATH_PRIM_META_SIZE_HPP
#define STAN_MATH_PRIM_META_SIZE_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
namespace math {
/** \ingroup type_trait
 * Returns the length of primitive scalar types
 * that are always of length 1.
 */
template <typename T, typename = require_stan_scalar_t<T>>
inline size_t size(const T& /*x*/) {
  return 1U;
}

/** \ingroup type_trait
 * Returns the size of the provided Eigen matrix, expression or std::vector.
 *
 * @param m input  \c Eigen \c Matrix, expression or std::vector
 * @tparam T type of m
 */
template <typename T, typename = require_not_stan_scalar_t<T>, typename = void>
inline size_t size(const T& m) {
  return m.size();
}
}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_META_LENGTH_HPP
#define STAN_MATH_PRIM_META_LENGTH_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
/** \ingroup type_trait
 * Returns the length of primitive scalar types
 * that are always of length 1.
 */
template <typename T, typename = require_stan_scalar_t<T>>
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

/** \ingroup type_trait
 * Returns the size of the provided Eigen matrix.
 *
 * @param m a const Eigen matrix
 * @tparam T type of matrix.
 * @return the size of the input matrix
 */
template <typename T, typename = require_eigen_t<T>, typename = void>
size_t length(const T& m) {
  return m.size();
}
}  // namespace stan
#endif

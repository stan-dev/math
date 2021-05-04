#ifndef STAN_MATH_PRIM_FUN_SIZE_HPP
#define STAN_MATH_PRIM_FUN_SIZE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Returns the length of primitive scalar types
 * that are always of length 1.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline size_t size(const T& /*x*/) {
  return 1U;
}

/** \ingroup type_trait
 * Returns the size of the provided Eigen matrix, expression or std::vector.
 *
 * @param m input  \c Eigen \c Matrix, expression or std::vector
 * @tparam T type of m
 */
template <typename T, require_container_t<T>* = nullptr>
inline size_t size(const T& m) {
  return m.size();
}

template <typename T, require_var_matrix_t<T>* = nullptr>
inline size_t size(const T& m) {
  return m.size();
}

}  // namespace math
}  // namespace stan
#endif

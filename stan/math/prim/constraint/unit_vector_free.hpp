#ifndef STAN_MATH_PRIM_CONSTRAINT_UNIT_VECTOR_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_UNIT_VECTOR_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Transformation of a unit length vector to a "free" vector
 * However, we are just fixing the unidentified radius to 1.
 * Thus, the transformation is just the identity
 *
 * @tparam EigVec A type derived from `EigenBase` with 1 compile time row or
 * column.
 * @param x unit vector of dimension K
 * @return Unit vector of dimension K considered "free"
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
inline auto unit_vector_free(EigVec&& x) {
  auto&& x_ref = to_ref(std::forward<EigVec>(x));
  check_unit_vector("stan::math::unit_vector_free", "Unit vector variable",
                    x_ref);
  return x_ref;
}

/**
 * Overload of `unit_vector_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto unit_vector_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return unit_vector_free(v); });
}

}  // namespace math
}  // namespace stan

#endif

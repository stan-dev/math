#ifndef STAN_MATH_PRIM_FUN_DETERMINANT_HPP
#define STAN_MATH_PRIM_FUN_DETERMINANT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the determinant of the specified square matrix.
 *
 * @tparam T type of the matrix (must be derived from \c Eigen::MatrixBase)
 *
 * @param m Specified matrix.
 * @return Determinant of the matrix.
 * @throw std::domain_error if matrix is not square.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline value_type_t<T> determinant(const T& m) {
  check_square("determinant", "m", m);
  return m.determinant();
}

}  // namespace math
}  // namespace stan

#endif

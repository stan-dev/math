#ifndef STAN_MATH_PRIM_MAT_FUN_DETERMINANT_HPP
#define STAN_MATH_PRIM_MAT_FUN_DETERMINANT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the determinant of the specified square matrix.
 *
 * @param m Specified matrix.
 * @return Determinant of the matrix.
 * @throw std::domain_error if matrix is not square.
 */
template <typename T, enable_if_eigen<T>* = nullptr>
inline auto determinant(const T& m) {
  check_square("determinant", "m", m);
  return m.determinant();
}

}  // namespace math
}  // namespace stan
#endif

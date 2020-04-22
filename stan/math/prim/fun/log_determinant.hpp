#ifndef STAN_MATH_PRIM_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_PRIM_FUN_LOG_DETERMINANT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the log absolute determinant of the specified square matrix.
 *
 * @tparam EigMat type of the matrix
 *
 * @param m Specified matrix.
 * @return log absolute determinant of the matrix.
 * @throw std::domain_error if matrix is not square.
 */
template <typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline value_type_t<EigMat> log_determinant(const EigMat& m) {
  if (m.size() == 0) {
    return 0;
  }
  check_square("log_determinant", "m", m);
  return m.colPivHouseholderQr().logAbsDeterminant();
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_PRIM_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_PRIM_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log absolute determinant of the specified square matrix.
 *
 * @tparam T type of the matrix
 *
 * @param m specified matrix
 * @return log absolute determinant of the matrix
 * @throw std::domain_error if matrix is not square and symmetric
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline value_type_t<EigMat> log_determinant_spd(const EigMat& m) {
  const auto& m_ref = to_ref(m);
  check_symmetric("log_determinant_spd", "m", m_ref);
  if (m.size() == 0) {
    return 0;
  }
  return sum(log(m_ref.ldlt().vectorD().array()));
}

}  // namespace math
}  // namespace stan

#endif

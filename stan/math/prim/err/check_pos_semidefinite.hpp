#ifndef STAN_MATH_PRIM_ERR_CHECK_POS_SEMIDEFINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POS_SEMIDEFINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <sstream>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is positive definite
 * @tparam T_y scalar type of the matrix
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::invalid_argument</code> if the matrix is not square
 *   or if the matrix has 0 size.
 * @throw <code>std::domain_error</code> if the matrix is not symmetric,
 *   or if it is not positive semi-definite,
 *   or if any element of the matrix is <code>NaN</code>.
 */
template <typename T_y>
inline void check_pos_semidefinite(
    const char* function, const char* name,
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  check_symmetric(function, name, y);
  check_positive(function, name, "rows", y.rows());
  check_not_nan(function, name, y);

  if (y.rows() == 1 && !(y(0, 0) >= 0.0)) {
    throw_domain_error(function, name, "is not positive semi-definite.", "");
  }

  using Eigen::Dynamic;
  using Eigen::LDLT;
  using Eigen::Matrix;
  LDLT<Matrix<double, Dynamic, Dynamic> > cholesky = value_of_rec(y).ldlt();
  if (cholesky.info() != Eigen::Success
      || (cholesky.vectorD().array() < 0.0).any()) {
    throw_domain_error(function, name, "is not positive semi-definite.", "");
  }
}

/**
 * Check if the specified matrix is positive semidefinite
 *
 * @tparam Derived Derived type of the Eigen::LDLT transform.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param cholesky Eigen::LDLT to test
 * @throw <code>std::domain_error</code> if the matrix is not positive
 *   semi-definite.
 */
template <typename Derived>
inline void check_pos_semidefinite(const char* function, const char* name,
                                   const Eigen::LDLT<Derived>& cholesky) {
  if (cholesky.info() != Eigen::Success
      || (cholesky.vectorD().array() < 0.0).any()) {
    throw_domain_error(function, name, "is not positive semi-definite.", "");
  }
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_ERR_CHECK_POS_DEFINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POS_DEFINITE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Check if the specified square, symmetric matrix is positive definite.
 *
 * This computes an LDLT decomposition to establish positive definiteness,
 * so it can be expensive. If an LDLT or LLT decomposition is available,
 * that should be tested instead.
 *
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param y matrix to test
 * @throw std::invalid_argument if the matrix is not square
 * or if the matrix has 0 size.
 * @throw std::domain_error if the matrix is not symmetric,
 * if it is not positive definite, or if any element is NaN
 */
template <typename EigMat, require_matrix_t<EigMat>* = nullptr>
inline void check_pos_definite(const char* function, const char* name,
                               const EigMat& y) {
  const auto& y_ref = to_ref(value_of_rec(y));
  check_symmetric(function, name, y_ref);
  check_positive(function, name, "rows", y_ref.rows());
  check_not_nan(function, name, y_ref);

  if (y_ref.rows() == 1 && !(y_ref(0, 0) > CONSTRAINT_TOLERANCE)) {
    throw_domain_error(function, name, "is not positive definite.", "");
  }

  Eigen::LDLT<Eigen::MatrixXd> cholesky = value_of_rec(y_ref).ldlt();
  if (cholesky.info() != Eigen::Success || !cholesky.isPositive()
      || (cholesky.vectorD().array() <= 0.0).any()) {
    throw_domain_error(function, name, "is not positive definite.", "");
  }
}

/**
 * Check if the specified LDLT decomposition of a matrix is positive definite.
 *
 * @tparam Derived type of the Eigen::LDLT decomposition
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param cholesky Eigen::LDLT to test, whose progenitor
 * must not have any NaN elements
 * @throw std::domain_error if the matrix is not positive definite
 */
template <typename Derived>
inline void check_pos_definite(const char* function, const char* name,
                               const Eigen::LDLT<Derived>& cholesky) {
  if (cholesky.info() != Eigen::Success || !cholesky.isPositive()
      || !(cholesky.vectorD().array() > 0.0).all()) {
    throw_domain_error(function, "LDLT decomposition of", " failed", name);
  }
}

/**
 * Check if the specified LLT decomposition was successful.
 *
 * @tparam Derived type of the Eigen::LLT decomposition
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param cholesky Eigen::LLT to test, whose progenitor
 * must not have any NaN elements
 * @throw std::domain_error if the decomposition failed or the
 * diagonal of the L matrix is not positive
 */
template <typename Derived>
inline void check_pos_definite(const char* function, const char* name,
                               const Eigen::LLT<Derived>& cholesky) {
  if (cholesky.info() != Eigen::Success
      || !(cholesky.matrixLLT().diagonal().array() > 0.0).all()) {
    throw_domain_error(function, "Matrix", " is not positive definite", name);
  }
}

}  // namespace math
}  // namespace stan
#endif

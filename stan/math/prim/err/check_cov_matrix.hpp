#ifndef STAN_MATH_PRIM_ERR_CHECK_COV_MATRIX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_COV_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/err/check_pos_definite.hpp>

namespace stan {
namespace math {
/**
 * Throw an exception if the specified matrix is not a valid covariance matrix.
 * A valid covariance matrix is a square, symmetric matrix that is positive
 * definite.
 * @tparam Mat Type inheriting from `MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1 or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw `std::invalid_argument` if the matrix is not square or if the matrix
 * is 0x0
 * @throw `std::domain_error` if the matrix is not symmetric, if the matrix is
 * not positive definite, or if any element of the matrix is `NaN`
 */
template <typename Mat, require_matrix_t<Mat>* = nullptr>
inline void check_cov_matrix(const char* function, const char* name,
                             const Mat& y) {
  check_pos_definite(function, name, y);
}

/**
 * Throw an exception if the specified matrix is not a valid covariance matrix.
 * A valid covariance matrix is a square, symmetric matrix that is positive
 * definite.
 * @tparam StdVec A standard vector with inner type either inheriting from
 * `MatrixBase` with neither rows or columns defined at compile time to be equal
 * to 1 or a `var_value` with the var's inner type inheriting from
 * `Eigen::MatrixBase` with neither rows or columns defined at compile time to
 * be equal to 1
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y standard vector of matrices to test.
 * @throw `std::invalid_argument` if the matrix is not square
 *   or if the matrix is 0x0
 * @throw `std::domain_error` if the matrix is not symmetric, if the matrix is
 * not positive definite, or if any element of the matrix is `NaN`
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
void check_cov_matrix(const char* function, const char* name, const StdVec& y) {
  for (auto&& y_i : y) {
    check_cov_matrix(function, name, y_i);
  }
}

}  // namespace math
}  // namespace stan
#endif

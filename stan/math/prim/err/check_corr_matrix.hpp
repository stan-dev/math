#ifndef STAN_MATH_PRIM_ERR_CHECK_CORR_MATRIX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CORR_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/check_pos_definite.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <sstream>
#include <string>
#include <cmath>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified matrix is not a valid correlation matrix.
 * A valid correlation matrix is symmetric positive definite, has a unit
 * diagonal (all 1 values), and has all values between -1 and 1 (inclusive).
 * @tparam Mat Type inheriting from `MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1 or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1
 * @param function Name of the function this was called from
 * @param name Name of the variable
 * @param y Matrix to test
 * @throw `std::invalid_argument` if the matrix is not square
 * @throw `std::domain_error` if the matrix is non-symmetric, diagonals not near
 * 1, not positive definite, or any of the elements are `NaN`
 */
template <typename Mat, require_matrix_t<Mat>* = nullptr>
inline void check_corr_matrix(const char* function, const char* name,
                              const Mat& y) {
  auto&& y_ref = to_ref(value_of_rec(y));
  check_square(function, name, y_ref);
  using std::fabs;
  if (y_ref.size() == 0) {
    return;
  }

  for (Eigen::Index k = 0; k < y.rows(); ++k) {
    if (!(fabs(y_ref.coeff(k, k) - 1.0) <= CONSTRAINT_TOLERANCE)) {
      [&y_ref, name, k, function]() STAN_COLD_PATH {
        std::ostringstream msg;
        msg << "is not a valid correlation matrix. " << name << "("
            << stan::error_index::value + k << ","
            << stan::error_index::value + k << ") is ";
        std::string msg_str(msg.str());
        throw_domain_error(function, name, y_ref(k, k), msg_str.c_str(),
                           ", but should be near 1.0");
      }();
    }
  }
  check_pos_definite(function, name, y_ref);
}

/**
 * Throw an exception if the specified matrix is not a valid correlation matrix.
 * A valid correlation matrix is symmetric positive definite, has a unit
 * diagonal (all 1 values), and has all values between -1 and 1 (inclusive).
 * @tparam StdVec A standard vector with inner type either inheriting from
 * `Eigen::MatrixBase` with neither rows or columns defined at compile time to
 * be equal to 1 or a `var_value` with the var's inner type inheriting from
 * `Eigen::MatrixBase` with neither rows or columns defined at compile time to
 * be equal to 1.
 * @param function Name of the function this was called from
 * @param name Name of the variable
 * @param y Matrix to test
 * @throw `std::invalid_argument` if the matrix is not square
 * @throw `std::domain_error` if the matrix is non-symmetric, diagonals not near
 * 1, not positive definite, or any of the elements are `NaN`
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
void check_corr_matrix(const char* function, const char* name,
                       const StdVec& y) {
  for (auto&& y_i : y) {
    check_corr_matrix(function, name, y_i);
  }
}

}  // namespace math
}  // namespace stan
#endif

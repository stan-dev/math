#ifndef STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_CORR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_CORR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/err/check_lower_triangular.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/err/check_unit_vector.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified matrix is not a valid Cholesky factor of
 * a correlation matrix. A Cholesky factor is a lower triangular matrix whose
 * diagonal elements are all positive and each row has unit Euclidean length.
 * Note that Cholesky factors need not be square, but require at least as many
 * rows M as columns N (i.e., `M >= N`). Tolerance is specified by
 * `math::CONSTRAINT_TOLERANCE`. Tolerance is specified by
 * `math::CONSTRAINT_TOLERANCE`.
 * @tparam Mat Type inheriting from `MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1 or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase` with neither rows or columns
 * defined at compile time to be equal to 1
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw `std::domain_error` if y is not a valid Cholesky factor, if number of
 * rows is less than the number of columns, if there are 0 columns, or if any
 * element in matrix is NaN
 */
template <typename Mat, require_matrix_t<Mat>* = nullptr>
void check_cholesky_factor_corr(const char* function, const char* name,
                                const Mat& y) {
  const auto& y_ref = to_ref(value_of_rec(y));
  check_square(function, name, y_ref);
  check_lower_triangular(function, name, y_ref);
  check_positive(function, name, y_ref.diagonal());
  for (Eigen::Index i = 0; i < y_ref.rows(); ++i) {
    check_unit_vector(function, name, y_ref.row(i));
  }
}

/**
 * Throw an exception if the specified matrix is not a valid Cholesky factor of
 * a correlation matrix. A Cholesky factor is a lower triangular matrix whose
 * diagonal elements are all positive and each row has unit Euclidean length.
 * Note that Cholesky factors need not be square, but require at least as many
 * rows M as columns N (i.e., `M >= N`). Tolerance is specified by
 * `math::CONSTRAINT_TOLERANCE`. Tolerance is specified by
 * `math::CONSTRAINT_TOLERANCE`.
 * @tparam StdVec A standard vector with inner type either inheriting from
 * `MatrixBase` with neither rows or columns defined at compile time to be equal
 * to 1 or a `var_value` with the var's inner type inheriting from
 * `Eigen::MatrixBase` with neither rows or columns defined at compile time to
 * be equal to 1
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Standard vector of matrics to test
 * @throw `std::domain_error` if y[i] is not a valid Cholesky factor, if number
 * of rows is less than the number of columns, if there are 0 columns, or if any
 * element in matrix is NaN
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
void check_cholesky_factor_corr(const char* function, const char* name,
                                const StdVec& y) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_cholesky_factor_corr(function,
                               internal::make_iter_name(name, i).c_str(), y[i]);
  }
}
}  // namespace math
}  // namespace stan
#endif

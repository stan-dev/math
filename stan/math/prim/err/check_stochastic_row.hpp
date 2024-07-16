#ifndef STAN_MATH_PRIM_ERR_CHECK_STOCHASTIC_ROW_HPP
#define STAN_MATH_PRIM_ERR_CHECK_STOCHASTIC_ROW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified matrix is not a row stochastic matrix. To
 * be a row stochastic matrix, all the values in each row must be greater than
 * or equal to 0 and the values must sum to 1. A valid row stochastic matrix is
 * one where the sum of the elements by row is equal to 1.  This function tests
 * that the sum is within the tolerance specified by `CONSTRAINT_TOLERANCE`.
 * This function only accepts Eigen matrices, statically typed vectors, not
 * general matrices with 1 column.
 * @tparam T A type inheriting from `Eigen::EigenBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Matrix to test
 * @throw `std::invalid_argument` if `theta` is a 0-vector
 * @throw `std::domain_error` if the vector is not a row stochastic matrix or if
 * any element is `NaN`
 */
template <typename T, require_matrix_t<T>* = nullptr>
void check_stochastic_row(const char* function, const char* name,
                          const T& theta) {
  using std::fabs;
  check_nonzero_size(function, name, theta);
  auto&& theta_ref = to_ref(value_of_rec(theta));
  for (Eigen::Index i = 0; i < theta_ref.rows(); ++i) {
    value_type_t<decltype(theta_ref)> vec_sum = 0.0;
    for (Eigen::Index j = 0; j < theta_ref.cols(); ++j) {
      if (!(theta_ref.coeff(i, j) >= 0)) {
        [&]() STAN_COLD_PATH {
          std::ostringstream msg;
          msg << "is not a valid row stochastic matrix. " << name << "["
              << std::to_string(i + stan::error_index::value) << ", "
              << std::to_string(i + stan::error_index::value) << "]"
              << " = ";
          std::string msg_str(msg.str());
          throw_domain_error(function, name, theta_ref.coeff(i, j),
                             msg_str.c_str(),
                             ", but should be greater than or equal to 0");
        }();
      }
      vec_sum += theta_ref.coeff(i, j);
    }
    if (!(fabs(1.0 - vec_sum) <= CONSTRAINT_TOLERANCE)) {
      [&]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << "is not a valid row stochastic matrix.";
        msg.precision(10);
        msg << " sum(" << name << "[" << std::to_string(i + 1)
            << ",:]) = " << vec_sum << ", but should be ";
        std::string msg_str(msg.str());
        throw_domain_error(function, name, 1.0, msg_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if the specified matrices in a standard vector are not a
 * row stochastic matrix. To be a row stochastic matrix, all the values in each
 * row must be greater than or equal to 0 and the values must sum to 1. A valid
 * row stochastic matrix is one where the sum of the elements by row is equal
 * to 1.  This function tests that the sum is within the tolerance specified by
 * `CONSTRAINT_TOLERANCE`. This function only accepts Eigen matrices, statically
 * typed vectors, not general matrices with 1 column.
 * @tparam T A type inheriting from `Eigen::EigenBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Matrix to test
 * @throw `std::invalid_argument` if `theta` is a 0-vector
 * @throw `std::domain_error` if the standard vector's matrices are not row
 * stochastic matrix or if any element is `NaN`
 */
template <typename T, require_std_vector_t<T>* = nullptr>
void check_stochastic_row(const char* function, const char* name,
                          const T& theta) {
  for (size_t i = 0; i < theta.size(); ++i) {
    check_stochastic_row(function, internal::make_iter_name(name, i).c_str(),
                         theta[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

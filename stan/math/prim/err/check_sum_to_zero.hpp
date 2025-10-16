#ifndef STAN_MATH_PRIM_ERR_CHECK_SUM_TO_ZERO_HPP
#define STAN_MATH_PRIM_ERR_CHECK_SUM_TO_ZERO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
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
 * Throw an exception if the specified vector does not sum to 0.
 * This function tests that the sum is within the tolerance specified by
 * `CONSTRAINT_TOLERANCE`.
 * This function only accepts Eigen vectors, statically
 * typed vectors, not general matrices with 1 column.
 * @tparam T A type inheriting from `Eigen::EigenBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test
 * @throw `std::invalid_argument` if `theta` is a 0-vector
 * @throw `std::domain_error` if the vector does not sum to zero
 */
template <typename T, require_matrix_t<T>* = nullptr>
void check_sum_to_zero(const char* function, const char* name, const T& theta) {
  using std::fabs;
  // the size-zero case is technically a valid sum-to-zero vector,
  // but it cannot be unconstrained to anything
  check_nonzero_size(function, name, theta);
  auto&& theta_ref = to_ref(value_of_rec(theta));
  if (unlikely(!(fabs(theta_ref.sum()) <= CONSTRAINT_TOLERANCE))) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      scalar_type_t<T> sum = theta_ref.sum();
      msg << "does not sum to zero.";
      msg.precision(10);
      msg << " sum(" << name << ") = " << sum << ", but should be ";
      std::string msg_str(msg.str());
      throw_domain_error(function, name, 0.0, msg_str.c_str());
    }();
  }
}

/**
 * Throw an exception if any vector in a standard vector does not sum to 0.
 * This function tests that the sum is within the tolerance specified by
 * `CONSTRAINT_TOLERANCE`.
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::EigenBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test.
 * @throw `std::invalid_argument` if `theta` is a 0-vector
 * @throw `std::domain_error` if the vector does not sum to zero
 */
template <typename T, require_std_vector_t<T>* = nullptr>
void check_sum_to_zero(const char* function, const char* name, const T& theta) {
  for (size_t i = 0; i < theta.size(); ++i) {
    check_sum_to_zero(function, internal::make_iter_name(name, i).c_str(),
                      theta[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

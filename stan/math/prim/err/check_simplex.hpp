#ifndef STAN_MATH_PRIM_ERR_CHECK_SIMPLEX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_SIMPLEX_HPP

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
 * Check if the specified vector is simplex.
 * To be a simplex, all values must be greater than or equal to 0
 * and the values must sum to 1.
 * A valid simplex is one where the sum of the elements is equal
 * to 1.  This function tests that the sum is within the tolerance
 * specified by <code>CONSTRAINT_TOLERANCE</code>. This function
 * only accepts Eigen vectors, statically typed vectors, not
 * general matrices with 1 column.
 * @tparam T A type inheriting from EigenBase.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test.
 * @throw <code>std::invalid_argument</code> if <code>theta</code>
 *   is a 0-vector.
 * @throw <code>std::domain_error</code> if the vector is not a
 *   simplex or if any element is <code>NaN</code>.
 */
template <typename T, require_matrix_t<T>* = nullptr>
void check_simplex(const char* function, const char* name, const T& theta) {
  using std::fabs;
  check_nonzero_size(function, name, theta);
  auto&& theta_ref = to_ref(value_of_rec(theta));
  if (!(fabs(1.0 - theta_ref.sum()) <= CONSTRAINT_TOLERANCE)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      value_type_t<T> sum = theta_ref.sum();
      msg << "is not a valid simplex.";
      msg.precision(10);
      msg << " sum(" << name << ") = " << sum << ", but should be ";
      std::string msg_str(msg.str());
      throw_domain_error(function, name, 1.0, msg_str.c_str());
    }();
  }
  for (Eigen::Index n = 0; n < theta_ref.size(); n++) {
    if (!(theta_ref.coeff(n) >= 0)) {
      [&]() STAN_COLD_PATH {
        std::ostringstream msg;
        msg << "is not a valid simplex. " << name << "["
            << n + stan::error_index::value << "]"
            << " = ";
        std::string msg_str(msg.str());
        throw_domain_error(function, name, theta_ref.coeff(n), msg_str.c_str(),
                           ", but should be greater than or equal to 0");
      }();
    }
  }
}

/**
 * Check if each vector in a standard vector is a simplex.
 * @tparam T A standard vector with inner type inheriting from `EigenBase`.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test.
 * @throw <code>std::invalid_argument</code> if <code>theta</code>
 *   is a 0-vector.
 * @throw <code>std::domain_error</code> if the vector is not a
 *   simplex or if any element is <code>NaN</code>.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
void check_simplex(const char* function, const char* name, const T& theta) {
  for (size_t i = 0; i < theta.size(); ++i) {
    check_simplex(function, internal::make_iter_name(name, i).c_str(),
                  theta[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

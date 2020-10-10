#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_scal_finite.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {
/**
 * Return true if y is finite
 *
 * @tparam T_y type of y
 * @param y parameter to check
 * @return boolean
 */
template <typename T_y>
bool is_finite(const T_y& y) {
  return is_scal_finite(y);
}

/**
 * Return true if every element of the matrix y is finite
 *
 * @tparam T_y type of elements y
 * @param y matrix to check
 * @return boolean
 */
template <typename T_y, int R, int C>
bool is_finite(const Eigen::Matrix<T_y, R, C>& y) {
  bool all = true;
  for (size_t n = 0; n < y.size(); ++n) {
    all &= is_finite(y(n));
  }
  return all;
}

/**
 * Return true if every element of the vector y is finite
 *
 * @tparam T_y type of elements y
 * @param y vector to check
 * @return boolean
 */
template <typename T_y>
bool is_finite(const std::vector<T_y>& y) {
  bool all = true;
  for (size_t n = 0; n < stan::math::size(y); ++n) {
    all &= is_finite(y[n]);
  }
  return all;
}
}  // namespace internal

/**
 * Check if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if y is infinity, -infinity, or NaN
 */
template <typename T_y, require_stan_scalar_t<T_y>* = nullptr>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  if (!internal::is_finite(y)) {
    throw_domain_error(function, name, y, "is ", ", but must be finite!");
  }
}

/**
 * Return <code>true</code> if all values in the std::vector are finite.
 *
 * @tparam T_y type of elements in the std::vector
 *
 * @param function name of function (for error messages)
 * @param name variable name (for error messages)
 * @param y std::vector to test
 * @return <code>true</code> if all values are finite
 **/
template <typename T_y, require_stan_scalar_t<T_y>* = nullptr>
inline void check_finite(const char* function, const char* name,
                         const std::vector<T_y>& y) {
  for (size_t n = 0; n < stan::math::size(y); n++) {
    if (!internal::is_finite(stan::get(y, n))) {
      throw_domain_error_vec(function, name, y, n, "is ",
                             ", but must be finite!");
    }
  }
}

/**
 * Return <code>true</code> is the specified matrix is finite.
 *
 * @tparam Derived Eigen derived type
 *
 * @param function name of function (for error messages)
 * @param name variable name (for error messages)
 * @param y matrix to test
 * @return <code>true</code> if the matrix is finite
 **/
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline void check_finite(const char* function, const char* name,
                         const EigMat& y) {
  if (!value_of(y).allFinite()) {
    for (int n = 0; n < y.size(); ++n) {
      if (!std::isfinite(value_of_rec(y(n)))) {
        throw_domain_error_vec(function, name, y, n, "is ",
                               ", but must be finite!");
      }
    }
  }
}

/**
 * Return <code>true</code> if all values in the std::vector are finite.
 *
 * @tparam T_y type of elements in the std::vector
 *
 * @param function name of function (for error messages)
 * @param name variable name (for error messages)
 * @param y std::vector to test
 * @return <code>true</code> if all values are finite
 **/
template <typename T_y, require_not_stan_scalar_t<T_y>* = nullptr>
inline void check_finite(const char* function, const char* name,
                         const std::vector<T_y>& y) {
  for (size_t n = 0; n < stan::math::size(y); n++) {
    if (!internal::is_finite(stan::get(y, n))) {
      throw_domain_error(function, name, "", "", "is not finite!");
    }
  }
}

}  // namespace math
}  // namespace stan

#endif

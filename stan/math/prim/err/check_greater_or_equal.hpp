#ifndef STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is greater than or equal to <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y Type of y
 * @tparam T_low Type of Lower bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not less than low
 *   or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_all_stan_scalar_t<T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name, const T_y& y,
                       const T_low& low) {
  if (!(y >= low)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << ", but must be greater than or equal to";
      msg << low;
      std::string msg_str(msg.str());
      throw_domain_error(function, name, y, "is ", msg_str.c_str());
    }();
  }
}

/**
 * Check if <code>y</code> is greater than or equal to <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y A container type
 * @tparam T_low A container type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not less than low
 *   or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_all_container_t<T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name, const T_y& y,
                       const T_low& low) {
  const auto& low_ref = as_array_or_scalar(low);
  const auto& y_ref = as_array_or_scalar(y);
  Eigen::Index n = 0;
  for (Eigen::Index j = 0; j < y_ref.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_ref.rows(); ++i) {
      if (!(y_ref.coeff(i, j) >= low_ref.coeff(i, j))) {
        [&]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to";
          msg << low_ref.coeff(i, j);
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_ref.coeff(i, j), i + j,
                                 "is ", msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Check if <code>y</code> is greater than or equal to <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y A container type
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not less than low
 *   or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_container_t<T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name, const T_y& y,
                       const T_low& low) {
  const auto& y_ref = as_array_or_scalar(y);
  Eigen::Index n = 0;
  for (Eigen::Index j = 0; j < y_ref.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_ref.rows(); ++i) {
      if (!(y_ref.coeff(i, j) >= low)) {
        [&]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to";
          msg << low;
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_ref.coeff(i, j), i + j,
                                 "is ", msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Check if <code>y</code> is greater than or equal to <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y A container type
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not less than low
 *   or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_container_t<T_low>* = nullptr,
          require_stan_scalar_t<T_y>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name, const T_y& y,
                       const T_low& low) {
  const auto& low_ref = as_array_or_scalar(low);
  Eigen::Index n = 0;
  for (Eigen::Index j = 0; j < low_ref.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_ref.rows(); ++i) {
      if (!(y < low_ref.coeff(i, j))) {
        [&]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << low_ref(i, j);
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
  }
}

}  // namespace math
}  // namespace stan
#endif

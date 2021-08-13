#ifndef STAN_MATH_PRIM_ERR_CHECK_GREATER_HPP
#define STAN_MATH_PRIM_ERR_CHECK_GREATER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/err/throw_domain_error_mat.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Throw an exception if `y` is strictly greater than `low`. This function is vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y A scalar type type
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_all_stan_scalar_t<T_y, T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  check_not_nan(function, name, y);
  check_not_nan(function, "lower", low);
  if (!(y > low)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << ", but must be greater than ";
      msg << low;
      std::string msg_str(msg.str());
      throw_domain_error(function, name, y, "is ", msg_str.c_str());
    }();
  }
}

/**
 * Throw an exception if `y` is strictly greater than each element of `low`. This function is
 * vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y A scalar type
 * @tparam T_low Type inheriting from `MatrixBase` or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase`.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_stan_scalar_t<T_y>* = nullptr,
          require_matrix_t<T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  auto&& low_arr = as_array_or_scalar(to_ref(value_of_rec(low)));
  check_not_nan(function, name, y);
  check_not_nan(function, "lower", low_arr);
  for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
      if (!(y > low_arr.coeff(i, j))) {
        [&low_arr, y, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than ";
          msg << low_arr.coeff(i, j);
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is strictly greater than `low`. This function is
 * vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y Type inheriting from `MatrixBase` or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase`.
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_matrix_t<T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  auto&& y_arr = as_array_or_scalar(to_ref(value_of_rec(y)));
  check_not_nan(function, name, y_arr);
  check_not_nan(function, "lower", low);
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (!(y_arr.coeff(i, j) > low)) {
        [&y_arr, low, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than ";
          msg << low;
          std::string msg_str(msg.str());
          throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is strictly greater than the associated element in `low`. This function is vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y Type inheriting from `MatrixBase` or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase`.
 * @tparam T_low Type inheriting from `MatrixBase` or a `var_value` with the var's
 * inner type inheriting from `Eigen::MatrixBase`.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_all_matrix_t<T_y, T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  auto&& y_arr = to_ref(as_array_or_scalar(value_of_rec(y)));
  auto&& low_arr = to_ref(as_array_or_scalar(value_of_rec(low)));
  check_not_nan(function, name, y_arr);
  check_not_nan(function, "lower", low_arr);
  if (is_vector<T_y>::value && is_vector<T_low>::value) {
    check_matching_sizes(function, name, y_arr, "lower", low_arr);
    for (Eigen::Index i = 0; i < low_arr.size(); ++i) {
      if (!(y_arr.coeff(i) > low_arr.coeff(i))) {
        [&y_arr, &low_arr, name, function, i]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than ";
          msg << low_arr.coeff(i);
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_arr, i, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  } else {
    check_matching_dims(function, name, y_arr, "lower", low_arr);
    for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
      for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
        if (!(y_arr.coeff(i, j) > low_arr.coeff(i, j))) {
          [&y_arr, &low_arr, name, function, i, j]() STAN_COLD_PATH {
            std::stringstream msg;
            msg << ", but must be greater than ";
            msg << low_arr.coeff(i, j);
            std::string msg_str(msg.str());
            throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                   msg_str.c_str());
          }();
        }
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is strictly greater than `low`. This function is vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam T_low A standard vector type whose `value_type` is a scalar, type inheriting from `Eigen::EigenBase`, or another standard vector
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_all_std_vector_t<T_y, T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  check_matching_sizes(function, name, y, "lower", low);
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater(function, internal::make_iter_name(name, i).c_str(), y[i],
                  low[i]);
  }
}

/**
 * Throw an exception if each element of `y` is strictly greater than `low`. This function is vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam T_low A scalar type or the same type as the underlying type in `T_y`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_std_vector_t<T_y>* = nullptr,
          require_not_std_vector_t<T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater(function, internal::make_iter_name(name, i).c_str(), y[i],
                  low);
  }
}

/**
 * Throw an exception if each element of `y` is strictly greater than `low`. This function is vectorized and will check each element of `y` against each element of `low`.
 * @tparam T_y A scalar type or the same type as the inner type of `T_low`
 * @tparam T_low A standard vector type whose `value_type` is a scalar, type inheriting from `Eigen::EigenBase`, or another standard vector
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any element of y or low is NaN.
 */
template <typename T_y, typename T_low,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_t<T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  for (size_t i = 0; i < low.size(); ++i) {
    check_greater(function, name, y, low[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

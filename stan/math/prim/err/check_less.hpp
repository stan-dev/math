#ifndef STAN_MATH_PRIM_ERR_CHECK_LESS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LESS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_mat.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Throw an exception if `y` is not strictly less than `high`. This function is vectorized and
 * will check each element of `y` against each element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_stan_scalar_t<T_y, T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  check_not_nan(function, name, y);
  check_not_nan(function, "higher", high);
  if (!(y < high)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << ", but must be less than ";
      msg << high;
      std::string msg_str(msg.str());
      throw_domain_error(function, name, y, "is ", msg_str.c_str());
    }();
  }
}

/**
 * Throw an exception if `y` is not strictly less than each element of `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high Type inheriting from `MatrixBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::MatrixBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
          require_matrix_t<T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  auto&& high_arr = as_array_or_scalar(to_ref(value_of_rec(high)));
  check_not_nan(function, name, y);
  check_not_nan(function, "higher", high_arr);
  for (Eigen::Index j = 0; j < high_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < high_arr.rows(); ++i) {
      if (!(y < high_arr.coeff(i, j))) {
        [&high_arr, y, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << high_arr.coeff(i, j);
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_high Type inheriting from `MatrixBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::MatrixBase`
 * @tparam T_high A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high, require_matrix_t<T_y>* = nullptr,
          require_stan_scalar_t<T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(y));
  check_not_nan(function, name, y_arr);
  check_not_nan(function, "higher", high);
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (!(y_arr.coeff(i, j) < high)) {
        [&y_arr, high, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << high;
          std::string msg_str(msg.str());
          throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than each element of `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y Type inheriting from `Eigen::MatrixBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::MatrixBase`
 * @tparam T_high Type inheriting from `Eigen::MatrixBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::MatrixBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_matrix_t<T_y, T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(y));
  auto&& high_arr = to_ref(value_of_rec(high));
  check_not_nan(function, name, y_arr);
  check_not_nan(function, "higher", high_arr);
  if (is_vector<T_y>::value && is_vector<T_high>::value) {
    check_matching_sizes(function, name, y_arr, "higher", high_arr);
    for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
      if (!(y_arr.coeff(i) < high_arr.coeff(i))) {
        [&y_arr, &high_arr, name, function, i]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << high_arr.coeff(i);
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_arr, i, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  } else {
    check_matching_dims(function, name, y_arr, "higher", high_arr);
    for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
      for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
        if (!(y_arr.coeff(i, j) < high_arr.coeff(i, j))) {
          [&y_arr, &high_arr, name, function, i, j]() STAN_COLD_PATH {
            std::stringstream msg;
            msg << ", but must be less than ";
            msg << high_arr.coeff(i, j);
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
 * Throw an exception if each element of `y` is not strictly less than each associated element of
 * `high`. This function is vectorized and will check each element of `y`
 * against each element of `high`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam T_high A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_std_vector_t<T_y, T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  check_matching_sizes(function, name, y, "higher", high);
  for (size_t i = 0; i < y.size(); ++i) {
    check_less(function, internal::make_iter_name(name, i).c_str(), y[i],
               high[i]);
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam T_high A scalar type or the same type as the inner type of `T_high`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high, require_std_vector_t<T_y>* = nullptr,
          require_not_std_vector_t<T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_less(function, internal::make_iter_name(name, i).c_str(), y[i], high);
  }
}

/**
 * Throw an exception if `y` is not strictly less than each element of `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_y A scalar type or the same type as the inner type of `T_high`
 * @tparam T_high A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_t<T_high>* = nullptr>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high) {
  for (size_t i = 0; i < high.size(); ++i) {
    check_less(function, name, y, high[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

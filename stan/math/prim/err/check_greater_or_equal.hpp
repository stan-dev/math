#ifndef STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/err/throw_domain_error_mat.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Throw an exception if `y` is not greater or equal than `low`. This function
 * is vectorized and will check each element of `y` against each element of
 * `low`.
 * @tparam T_y A scalar type
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_stan_scalar_t<T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  check_not_nan(function, name, y);
  check_not_nan(function, "lower", low);
  if (!(y >= low)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << ", but must be greater than or equal to ";
      msg << low;
      std::string msg_str(msg.str());
      throw_domain_error(function, name, y, "is ", msg_str.c_str());
    }();
  }
}

/**
 * Throw an exception if `y` is not greater or equal than each element of `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A scalar type
 * @tparam T_low A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_stan_scalar_t<T_y>* = nullptr,
          require_vector_vt<is_stan_scalar, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& low_arr = to_ref(value_of_rec(as_array_or_scalar(low)));
  for (Eigen::Index i = 0; i < low_arr.size(); ++i) {
    if (!(y >= low_arr.coeff(i))) {
      [&low_arr, y, name, function, i]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be greater than or equal to ";
        msg << low_arr.coeff(i);
        std::string msg_str(msg.str());
        throw_domain_error(function, name, y, "is ", msg_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if `y` is not greater or equal than each element of `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A scalar type
 * @tparam T_low Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_stan_scalar_t<T_y>* = nullptr,
          require_dense_dynamic_t<T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& low_arr = to_ref(value_of_rec(low));
  for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
      if (!(y >= low_arr.coeff(i, j))) {
        [&low_arr, y, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to ";
          msg << low_arr.coeff(i, j);
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_vector_vt<is_stan_scalar, T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& y_arr = to_ref(value_of_rec(as_array_or_scalar(y)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (!(y_arr.coeff(i) >= low)) {
      [&y_arr, low, name, function, i]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be greater than or equal to ";
        msg << low;
        std::string msg_str(msg.str());
        throw_domain_error_vec(function, name, y_arr, i, "is ",
                               msg_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_low A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_dense_dynamic_t<T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& y_arr = to_ref(value_of_rec(y));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (!(y_arr.coeff(i, j) >= low)) {
        [&y_arr, low, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to ";
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
 * Throw an exception if each element of `y` is not greater or equal than the
 * associated element in `low`. This function is vectorized and will check each
 * element of `y` against each element of `low`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar
 * @tparam T_low A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_vector_vt<is_stan_scalar, T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& y_arr = to_ref(value_of_rec(as_array_or_scalar(y)));
  auto&& low_arr = to_ref(value_of_rec(as_array_or_scalar(low)));
  check_matching_sizes(function, name, y_arr, "lower", low_arr);
  for (Eigen::Index i = 0; i < low_arr.size(); ++i) {
    if (!(y_arr.coeff(i) >= low_arr.coeff(i))) {
      [&y_arr, &low_arr, name, function, i]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be greater than or equal to ";
        msg << low_arr.coeff(i);
        std::string msg_str(msg.str());
        throw_domain_error_vec(function, name, y_arr, i, "is ",
                               msg_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than the
 * associated element in `low`. This function is vectorized and will check each
 * element of `y` against each element of `low`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_low Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_dense_dynamic_t<T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  auto&& y_arr = to_ref(value_of_rec(y));
  auto&& low_arr = to_ref(value_of_rec(low));
  check_matching_dims(function, name, y_arr, "lower", low_arr);
  for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
      if (!(y_arr.coeff(i, j) >= low_arr.coeff(i, j))) {
        [&y_arr, &low_arr, name, function, i, j]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to ";
          msg << low_arr.coeff(i, j);
          std::string msg_str(msg.str());
          throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector type
 * @tparam T_low A standard vector type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_any_std_vector_vt<is_container, T_y, T_low>* = nullptr,
          require_all_std_vector_t<T_y, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  check_matching_sizes(function, name, y, "lower", low);
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater_or_equal(function, internal::make_iter_name(name, i).c_str(),
                           y[i], low[i]);
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector type with a `value_type` of a standard vector
 * or type inheriting from `Eigen::DenseBase`
 * @tparam T_low A standard vector type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_std_vector_vt<is_container, T_y>* = nullptr,
          require_not_std_vector_t<T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater_or_equal(function, internal::make_iter_name(name, i).c_str(),
                           y[i], low);
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector
 * @tparam T_low  A standard vector type with a `value_type` of a standard
 * vector or type inheriting from `Eigen::DenseBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is NaN
 */
template <typename T_y, typename T_low,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_vt<is_container, T_low>* = nullptr>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  for (size_t i = 0; i < low.size(); ++i) {
    check_greater_or_equal(function, name, y, low[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif

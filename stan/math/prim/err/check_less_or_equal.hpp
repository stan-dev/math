#ifndef STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
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
 * Throw an exception if `y` is not less than `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_stan_scalar_t<T_y, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  if (unlikely(!(y <= high))) {
    [&]() STAN_COLD_PATH {
      throw_domain_error(function, name, y, "is ",
                         (", but must be less than or equal to "
                          + std::to_string(value_of_rec(high)))
                             .c_str());
    }();
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
          require_vector_vt<is_stan_scalar, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& high_arr = to_ref(value_of_rec(as_array_or_scalar(high)));
  for (Eigen::Index i = 0; i < high_arr.size(); ++i) {
    if (unlikely(!(y <= high_arr.coeff(i)))) {
      [&high_arr, y, name, function, i]() STAN_COLD_PATH {
        throw_domain_error(function, name, y, "is ",
                           (", but must be less than or equal to "
                            + std::to_string(high_arr.coeff(i)))
                               .c_str());
      }();
    }
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
          require_dense_dynamic_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& high_arr = to_ref(value_of_rec(high));
  for (Eigen::Index j = 0; j < high_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < high_arr.rows(); ++i) {
      if (unlikely(!(y <= high_arr.coeff(i, j)))) {
        [&high_arr, y, name, function, i, j]() STAN_COLD_PATH {
          throw_domain_error(function, name, y, "is ",
                             (", but must be less than or equal to "
                              + std::to_string(high_arr.coeff(i, j)))
                                 .c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam T_high A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_vector_vt<is_stan_scalar, T_y>* = nullptr,
          require_stan_scalar_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(as_array_or_scalar(y)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) <= high))) {
      [&y_arr, high, name, function, i]() STAN_COLD_PATH {
        throw_domain_error_vec(function, name, y_arr, i, "is ",
                               (", but must be less than or equal to "
                                + std::to_string(value_of_rec(high)))
                                   .c_str());
      }();
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_high A scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_dense_dynamic_t<T_y>* = nullptr,
          require_stan_scalar_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(y));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) <= high))) {
        [&y_arr, high, name, function, i, j]() STAN_COLD_PATH {
          throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                 (", but must be less than or equal to "
                                  + std::to_string(value_of_rec(high)))
                                     .c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than the associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam T_high A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_vector_vt<is_stan_scalar, T_y, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(as_array_or_scalar(y)));
  auto&& high_arr = to_ref(value_of_rec(as_array_or_scalar(high)));
  check_matching_sizes(function, name, y_arr, "higher", high_arr);
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) <= high_arr.coeff(i)))) {
      [&y_arr, &high_arr, name, function, i]() STAN_COLD_PATH {
        throw_domain_error_vec(function, name, y_arr, i, "is ",
                               (", but must be less than or equal to "
                                + std::to_string(high_arr.coeff(i)))
                                   .c_str());
      }();
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than the associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_high Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_dense_dynamic_t<T_y, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  auto&& y_arr = to_ref(value_of_rec(y));
  auto&& high_arr = to_ref(value_of_rec(high));
  check_matching_dims(function, name, y_arr, "higher", high_arr);
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) <= high_arr.coeff(i, j)))) {
        [&y_arr, &high_arr, name, function, i, j]() STAN_COLD_PATH {
          throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                 (", but must be less than or equal to "
                                  + std::to_string(high_arr.coeff(i, j)))
                                     .c_str());
        }();
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than each associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam T_high A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_any_std_vector_vt<is_container, T_y, T_high>* = nullptr,
          require_all_std_vector_t<T_y, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  check_matching_sizes(function, name, y, "higher", high);
  for (size_t i = 0; i < y.size(); ++i) {
    check_less_or_equal(function, name, y[i], high[i]);
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A standard vector type with a `value_type` of a standard vector
 * or type inheriting from `Eigen::DenseBase`
 * @tparam T_high A scalar or the same `value_type` of `T_y`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_std_vector_vt<is_container, T_y>* = nullptr,
          require_not_std_vector_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  std::string iter_name{name};
  for (size_t i = 0; i < y.size(); ++i) {
    check_less_or_equal(function, name, y[i], high);
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type or the same `value_type` of `T_high`
 * @tparam T_high A standard vector type with a `value_type` of a standard
 * vector or type inheriting from `Eigen::DenseBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename T_high,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_vt<is_container, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  for (size_t i = 0; i < high.size(); ++i) {
    check_less_or_equal(function, name, y, high[i]);
  }
}
}  // namespace math
}  // namespace stan
#endif

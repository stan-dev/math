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
 * Throw an exception if `y` is not strictly less than `high`. This function is
 * vectorized and will check each element of `y` against each element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_stan_scalar_t<T_y, T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  if (!(y < high)) {
    [](auto y, auto high, auto function, auto name,
       auto... idxs) STAN_COLD_PATH {
      throw_domain_error(
          function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
          (", but must be less than " + std::to_string(value_of_rec(high)))
              .c_str());
    }(y, high, function, name, idxs...);
  }
}

/**
 * Throw an exception if `y` is not strictly less than each element of `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <
    typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
    require_vector_t<T_high>* = nullptr,
    require_not_std_vector_vt<is_container_or_var_matrix, T_high>* = nullptr,
    typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& high_arr = as_array_or_scalar(value_of_rec(to_ref(high)));
  for (Eigen::Index i = 0; i < high_arr.size(); ++i) {
    if (unlikely(!(y < high_arr.coeff(i)))) {
      [](auto y, auto&& high_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error(
            function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
            (", but must be less than " + std::to_string(high_arr.coeff(i)))
                .c_str());
      }(y, high_arr, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if `y` is not strictly less than each element of `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type
 * @tparam T_high Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
          require_dense_dynamic_t<T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& high_arr = value_of_rec(to_ref(high));
  for (Eigen::Index j = 0; j < high_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < high_arr.rows(); ++i) {
      if (unlikely(!(y < high_arr.coeff(i, j)))) {
        [](auto y, auto&& high_arr, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error(function,
                             internal::make_iter_name(name, idxs...).c_str(), y,
                             "is ",
                             (", but must be less than "
                              + std::to_string(high_arr.coeff(i, j)))
                                 .c_str());
        }(y, high_arr, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam T_high A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high, require_vector_t<T_y>* = nullptr,
          require_not_std_vector_vt<is_container_or_var_matrix, T_y>* = nullptr,
          require_stan_scalar_t<T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& y_arr = value_of_rec(as_array_or_scalar(to_ref(y)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) < high))) {
      [](auto&& y_arr, auto high, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error_vec(
            function, internal::make_iter_name(name, idxs...).c_str(), y_arr, i,
            "is ",
            (", but must be less than " + std::to_string(value_of_rec(high)))
                .c_str());
      }(y_arr, high, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_high A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_dense_dynamic_t<T_y>* = nullptr,
          require_stan_scalar_t<T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) < high))) {
        [](auto&& y_arr, auto high, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be less than " + std::to_string(value_of_rec(high)))
                  .c_str());
        }(y_arr, high, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than each
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam T_y A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam T_high A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_vector_t<T_y, T_high>* = nullptr,
          require_all_not_std_vector_vt<is_container_or_var_matrix, T_y,
                                        T_high>* = nullptr,
          typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  auto&& high_arr = value_of_rec(to_ref(high));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) < high_arr.coeff(i)))) {
      [](auto&& y_arr, auto&& high_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error_vec(
            function, internal::make_iter_name(name, idxs...).c_str(), y_arr, i,
            "is ",
            (", but must be less than " + std::to_string(high_arr.coeff(i)))
                .c_str());
      }(y_arr, high_arr, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than each
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam T_y Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam T_high Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_all_dense_dynamic_t<T_y, T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  auto&& high_arr = value_of_rec(to_ref(high));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) < high_arr.coeff(i, j)))) {
        [](auto&& y_arr, auto&& high_arr, auto name, auto function, auto i,
           auto j, auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be less than "
               + std::to_string(high_arr.coeff(i, j)))
                  .c_str());
        }(y_arr, high_arr, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A standard vector type with a `value_type` of a standard vector
 * or type inheriting from `Eigen::DenseBase`
 * @tparam T_high A scalar type or the same type as the inner type of `T_high`
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_std_vector_vt<is_container_or_var_matrix, T_y>* = nullptr,
          require_not_std_vector_t<T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_less(function, name, y[i], high, idxs..., i);
  }
}

/**
 * Throw an exception if `y` is not strictly less than each element of `high`.
 * This function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type or the same type as the inner type of `T_high`
 * @tparam T_high A standard vector type with a `value_type` of a standard
 * vector or type inheriting from `Eigen::DenseBase`
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_vt<is_container_or_var_matrix, T_high>* = nullptr,
          typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  for (size_t i = 0; i < high.size(); ++i) {
    check_less(function, name, y, high[i], idxs..., i);
  }
}

/**
 * Throw an exception if each element of `y` is not strictly less than each
 * associated element of `high`. This function is vectorized and will check each
 * element of `y` against each element of `high`.
 * @tparam T_y A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::DenseBase`, or another standard vector
 * @tparam T_high A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::DenseBase`, or another standard vector
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `domain_error` if y is not less than high or if any element of y or
 * high is `NaN`
 */
template <typename T_y, typename T_high,
          require_any_std_vector_vt<is_container_or_var_matrix, T_y,
                                    T_high>* = nullptr,
          require_all_std_vector_t<T_y, T_high>* = nullptr, typename... Idxs>
inline void check_less(const char* function, const char* name, const T_y& y,
                       const T_high& high, Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_less(function, name, y[i], high[i], idxs..., i);
  }
}

}  // namespace math
}  // namespace stan
#endif

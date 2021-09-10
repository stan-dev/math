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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_stan_scalar_t<T_y, T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  if (unlikely(!(y >= low))) {
    [](auto y, auto low, auto function, auto name,
       auto... idxs) STAN_COLD_PATH {
      throw_domain_error(
          function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
          (", but must be greater than or equal to "
           + std::to_string(value_of_rec(low)))
              .c_str());
    }(y, low, function, name, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <
    typename T_y, typename T_low, require_stan_scalar_t<T_y>* = nullptr,
    require_vector_t<T_low>* = nullptr,
    require_not_std_vector_vt<is_container_or_var_matrix, T_low>* = nullptr,
    typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& low_arr = value_of_rec(as_array_or_scalar(to_ref(low)));
  for (Eigen::Index i = 0; i < low_arr.size(); ++i) {
    if (unlikely(!(y >= low_arr.coeff(i)))) {
      [](auto y, auto&& low_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error(
            function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
            (", but must be greater than or equal to "
             + std::to_string(low_arr.coeff(i)))
                .c_str());
      }(y, low_arr, name, function, i, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_stan_scalar_t<T_y>* = nullptr,
          require_dense_dynamic_t<T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& low_arr = value_of_rec(to_ref(low));
  for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
      if (unlikely(!(y >= low_arr.coeff(i, j)))) {
        [](auto y, auto&& low_arr, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error(function,
                             internal::make_iter_name(name, idxs...).c_str(), y,
                             "is ",
                             (", but must be greater than or equal to "
                              + std::to_string(low_arr.coeff(i, j)))
                                 .c_str());
        }(y, low_arr, name, function, i, j, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_vector_t<T_y>* = nullptr,
          require_not_std_vector_vt<is_container_or_var_matrix, T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& y_arr = value_of_rec(as_array_or_scalar(to_ref(y)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) >= low))) {
      [](auto&& y_arr, auto low, auto name, auto function, auto i, auto... idxs)
          STAN_COLD_PATH {
            throw_domain_error_vec(
                function, internal::make_iter_name(name, idxs...).c_str(),
                y_arr, i, "is ",
                (", but must be greater than or equal to "
                 + std::to_string(value_of_rec(low)))
                    .c_str());
          }(y_arr, low, name, function, i, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low, require_dense_dynamic_t<T_y>* = nullptr,
          require_stan_scalar_t<T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) >= low))) {
        [](auto&& y_arr, auto low, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be greater than or equal to "
               + std::to_string(value_of_rec(low)))
                  .c_str());
        }(y_arr, low, name, function, i, j, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_vector_t<T_y, T_low>* = nullptr,
          require_all_not_std_vector_vt<is_container_or_var_matrix, T_y,
                                        T_low>* = nullptr,
          typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& y_arr = value_of_rec(as_array_or_scalar(to_ref(y)));
  auto&& low_arr = value_of_rec(as_array_or_scalar(to_ref(low)));
  for (Eigen::Index i = 0; i < low_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) >= low_arr.coeff(i)))) {
      [](auto&& y_arr, auto&& low_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error_vec(function,
                               internal::make_iter_name(name, idxs...).c_str(),
                               y_arr, i, "is ",
                               (", but must be greater than or equal to "
                                + std::to_string(low_arr.coeff(i)))
                                   .c_str());
      }(y_arr, low_arr, name, function, i, idxs...);
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
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_all_dense_dynamic_t<T_y, T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  auto&& low_arr = value_of_rec(to_ref(low));
  for (Eigen::Index j = 0; j < low_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < low_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) >= low_arr.coeff(i, j)))) {
        [](auto&& y_arr, auto&& low_arr, auto name, auto function, auto i,
           auto j, auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be greater than or equal to "
               + std::to_string(low_arr.coeff(i, j)))
                  .c_str());
        }(y_arr, low_arr, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector type with a `value_type` of a standard vector
 * or type inheriting from `Eigen::DenseBase`
 * @tparam T_low A standard vector type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_std_vector_vt<is_container_or_var_matrix, T_y>* = nullptr,
          require_not_std_vector_t<T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater_or_equal(function, name, y[i], low, idxs..., i);
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector
 * @tparam T_low  A standard vector type with a `value_type` of a standard
 * @tparam Idxs A parameter pack of Integral types
 * vector or type inheriting from `Eigen::DenseBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is NaN
 */
template <typename T_y, typename T_low,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_vt<is_container_or_var_matrix, T_low>* = nullptr,
          typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  for (size_t i = 0; i < low.size(); ++i) {
    check_greater_or_equal(function, name, y, low[i], idxs..., i);
  }
}

/**
 * Throw an exception if each element of `y` is not greater or equal than `low`.
 * This function is vectorized and will check each element of `y` against each
 * element of `low`.
 * @tparam T_y A standard vector type
 * @tparam T_low A standard vector type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not greater or equal to low or if any
 * element of y or low is `NaN`
 */
template <typename T_y, typename T_low,
          require_any_std_vector_vt<is_container_or_var_matrix, T_y,
                                    T_low>* = nullptr,
          require_all_std_vector_t<T_y, T_low>* = nullptr, typename... Idxs>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low,
                                   Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_greater_or_equal(function, name, y[i], low[i], idxs..., i);
  }
}

}  // namespace math
}  // namespace stan
#endif

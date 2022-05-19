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
 * @tparam ScalarY A scalar type
 * @tparam ScalarHigh A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename ScalarY, typename ScalarHigh,
          require_all_stan_scalar_t<ScalarY, ScalarHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const ScalarY& y, const ScalarHigh& high,
                                Idxs... idxs) {
  if (unlikely(!(y <= high))) {
    [](auto y, auto high, auto function, auto name,
       auto... idxs) STAN_COLD_PATH {
      throw_domain_error(
          function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
          (", but must be less than or equal to "
           + std::to_string(value_of_rec(high)))
              .c_str());
    }(y, high, function, name, idxs...);
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam ScalarY A scalar type
 * @tparam VecHigh A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <
    typename ScalarY, typename VecHigh, require_stan_scalar_t<ScalarY>* = nullptr,
    require_vector_t<VecHigh>* = nullptr,
    require_not_std_vector_vt<is_container_or_var_matrix, VecHigh>* = nullptr,
    typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const ScalarY& y, const VecHigh& high,
                                Idxs... idxs) {
  auto&& high_arr = value_of_rec(as_array_or_scalar(to_ref(high)));
  for (Eigen::Index i = 0; i < high_arr.size(); ++i) {
    if (unlikely(!(y <= high_arr.coeff(i)))) {
      [](auto y, auto&& high_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error(
            function, internal::make_iter_name(name, idxs...).c_str(), y, "is ",
            (", but must be less than or equal to "
             + std::to_string(high_arr.coeff(i)))
                .c_str());
      }(y, high_arr, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam ScalarY A scalar type
 * @tparam DenseHigh Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename ScalarY, typename DenseHigh, require_stan_scalar_t<ScalarY>* = nullptr,
          require_dense_dynamic_t<DenseHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const ScalarY& y, const DenseHigh& high,
                                Idxs... idxs) {
  auto&& high_arr = value_of_rec(to_ref(high));
  for (Eigen::Index j = 0; j < high_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < high_arr.rows(); ++i) {
      if (unlikely(!(y <= high_arr.coeff(i, j)))) {
        [](auto y, auto&& high_arr, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error(function,
                             internal::make_iter_name(name, idxs...).c_str(), y,
                             "is ",
                             (", but must be less than or equal to "
                              + std::to_string(high_arr.coeff(i, j)))
                                 .c_str());
        }(y, high_arr, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam VecY A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam ScalarHigh A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename VecY, typename ScalarHigh, require_vector_t<VecY>* = nullptr,
          require_not_std_vector_vt<is_container_or_var_matrix, VecY>* = nullptr,
          require_stan_scalar_t<ScalarHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const VecY& y, const ScalarHigh& high,
                                Idxs... idxs) {
  auto&& y_arr = value_of_rec(as_array_or_scalar(to_ref(y)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) <= high))) {
      [](auto&& y_arr, auto high, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error_vec(function,
                               internal::make_iter_name(name, idxs...).c_str(),
                               y_arr, i, "is ",
                               (", but must be less than or equal to "
                                + std::to_string(value_of_rec(high)))
                                   .c_str());
      }(y_arr, high, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam DenseY Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam ScalarHigh A scalar type
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename DenseY, typename ScalarHigh,
          require_dense_dynamic_t<DenseY>* = nullptr,
          require_stan_scalar_t<ScalarHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const DenseY& y, const ScalarHigh& high,
                                Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) <= high))) {
        [](auto&& y_arr, auto high, auto name, auto function, auto i, auto j,
           auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be less than or equal to "
               + std::to_string(value_of_rec(high)))
                  .c_str());
        }(y_arr, high, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than the associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam DenseY A standard vector or type inheriting from `Eigen::DenseBase` with
 *  compile time rows or columns equal to one and `value_type` equal to a stan
 * scalar.
 * @tparam VecHigh A standard vector or type inheriting from `Eigen::DenseBase`
 * with compile time rows or columns equal to one and `value_type` equal to a
 * stan scalar.
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename DenseY, typename VecHigh,
          require_all_vector_t<DenseY, VecHigh>* = nullptr,
          require_all_not_std_vector_vt<is_container_or_var_matrix, DenseY,
                                        VecHigh>* = nullptr,
          typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const DenseY& y, const VecHigh& high,
                                Idxs... idxs) {
  auto&& y_arr = value_of_rec(as_array_or_scalar(to_ref(y)));
  auto&& high_arr = value_of_rec(as_array_or_scalar(to_ref(high)));
  for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
    if (unlikely(!(y_arr.coeff(i) <= high_arr.coeff(i)))) {
      [](auto&& y_arr, auto&& high_arr, auto name, auto function, auto i,
         auto... idxs) STAN_COLD_PATH {
        throw_domain_error_vec(function,
                               internal::make_iter_name(name, idxs...).c_str(),
                               y_arr, i, "is ",
                               (", but must be less than or equal to "
                                + std::to_string(high_arr.coeff(i)))
                                   .c_str());
      }(y_arr, high_arr, name, function, i, idxs...);
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than the associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam DenseY Type inheriting from `Eigen::DenseBase` or a `var_value` with the
 * var's inner type inheriting from `Eigen::DenseBase` where the compile time
 * number of rows or columns is not equal to one
 * @tparam DenseHigh Type inheriting from `Eigen::DenseBase` or a `var_value` with
 * the var's inner type inheriting from `Eigen::DenseBase` where the compile
 * time number of rows or columns is not equal to one
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename DenseY, typename DenseHigh,
          require_all_dense_dynamic_t<DenseY, DenseHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const DenseY& y, const DenseHigh& high,
                                Idxs... idxs) {
  auto&& y_arr = value_of_rec(to_ref(y));
  auto&& high_arr = value_of_rec(to_ref(high));
  for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
    for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
      if (unlikely(!(y_arr.coeff(i, j) <= high_arr.coeff(i, j)))) {
        [](auto&& y_arr, auto&& high_arr, auto name, auto function, auto i,
           auto j, auto... idxs) STAN_COLD_PATH {
          throw_domain_error_mat(
              function, internal::make_iter_name(name, idxs...).c_str(), y_arr,
              i, j, "is ",
              (", but must be less than or equal to "
               + std::to_string(high_arr.coeff(i, j)))
                  .c_str());
        }(y_arr, high_arr, name, function, i, j, idxs...);
      }
    }
  }
}

/**
 * Throw an exception if each element of `y` is not less than `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam ContainerY A standard vector type with a `value_type` of a standard vector
 * or type inheriting from `Eigen::DenseBase`
 * @tparam T_high A scalar or the same `value_type` of `ContainerY`
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename ContainerY, typename T_high,
          require_std_vector_vt<is_container_or_var_matrix, ContainerY>* = nullptr,
          require_not_std_vector_t<T_high>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const ContainerY& y, const T_high& high,
                                Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_less_or_equal(function, name, y[i], high, idxs..., i);
  }
}

/**
 * Throw an exception if `y` is not less than each element of `high`. This
 * function is vectorized and will check each element of `y` against each
 * element of `high`.
 * @tparam T_y A scalar type or the same `value_type` of `ContainerHigh`
 * @tparam ContainerHigh A standard vector type with a `value_type` of a standard
 * vector or type inheriting from `Eigen::DenseBase`
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename T_y, typename ContainerHigh,
          require_not_std_vector_t<T_y>* = nullptr,
          require_std_vector_vt<is_container_or_var_matrix, ContainerHigh>* = nullptr,
          typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const ContainerHigh& high,
                                Idxs... idxs) {
  for (size_t i = 0; i < high.size(); ++i) {
    check_less_or_equal(function, name, y, high[i], idxs..., i);
  }
}

/**
 * Throw an exception if each element of `y` is not less than each associated
 * element of `high`. This function is vectorized and will check each element of
 * `y` against each element of `high`.
 * @tparam StdVecY A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam StdVecHigh A standard vector type whose `value_type` is a scalar, type
 * inheriting from `Eigen::EigenBase`, or another standard vector
 * @tparam Idxs A parameter pack of Integral types
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @param idxs Pack of integral types to construct lazily construct the error
 * message indices
 * @throw `std::domain_error` if y is not less than high or if any element of y
 * or high is `NaN`
 */
template <typename StdVecY, typename StdVecHigh,
          require_any_std_vector_vt<is_container_or_var_matrix, StdVecY,
                                    StdVecHigh>* = nullptr,
          require_all_std_vector_t<StdVecY, StdVecHigh>* = nullptr, typename... Idxs>
inline void check_less_or_equal(const char* function, const char* name,
                                const StdVecY& y, const StdVecHigh& high,
                                Idxs... idxs) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_less_or_equal(function, name, y[i], high[i], idxs..., i);
  }
}

}  // namespace math
}  // namespace stan
#endif

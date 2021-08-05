#ifndef STAN_MATH_PRIM_ERR_CHECK_LESS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LESS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is strictly less than <code>high</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>high</code>.
 * @tparam T_y Type of y
 * @tparam T_high Type of upper bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw <code>domain_error</code> if y is not less than high
 *   or if any element of y or high is NaN.
 */
 template <typename T_y, typename T_high, require_all_stan_scalar_t<T_y, T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
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

 template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr, require_eigen_t<T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
    auto&& high_arr = as_array_or_scalar(to_ref(value_of_rec(high)));
    check_not_nan(function, name, y);
    check_not_nan(function, "higher", high_arr);
    for (Eigen::Index i = 0; i < high_arr.size(); ++i) {
      if (!(y < high_arr.coeffRef(i))) {
        [&high_arr, y, name, function, i]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << high_arr.coeff(i);
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
 }

 template <typename T_y, typename T_high, require_eigen_t<T_y>* = nullptr, require_stan_scalar_t<T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
    auto&& y_arr = as_array_or_scalar(to_ref(value_of_rec(y)));
    check_not_nan(function, name, y_arr);
    check_not_nan(function, "higher", high);
    for (Eigen::Index i = 0; i < y_arr.size(); ++i) {
      if (!(y_arr.coeffRef(i) < high)) {
        [&y_arr, high, name, function, i]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be less than ";
          msg << high;
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_arr, i, "is ",
                                 msg_str.c_str());
        }();
      }
    }
 }

 template <typename T_y, typename T_high, require_all_eigen_t<T_y, T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
    check_matching_dims(function, name, y, "higher", high);
    auto&& y_arr = as_array_or_scalar(to_ref(value_of_rec(y)));
    auto&& high_arr = as_array_or_scalar(to_ref(value_of_rec(high)));
    check_not_nan(function, name, y_arr);
    check_not_nan(function, "higher", high_arr);
    for (Eigen::Index i = 0; i < high_arr.size(); ++i) {
      if (!(y_arr.coeffRef(i) < high_arr.coeffRef(i))) {
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
 }

 template <typename T_y, typename T_high, require_all_std_vector_t<T_y, T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
   check_matching_dims(function, name, y, "higher", high);
   for (size_t i = 0; i < y.size(); ++i) {
     check_less(function, name, y[i], high[i]);
   }
 }

 template <typename T_y, typename T_high, require_std_vector_t<T_y>* = nullptr,
  require_not_std_vector_t<T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
   for (size_t i = 0; i < y.size(); ++i) {
     check_less(function, name, y[i], high);
   }
 }

 template <typename T_y, typename T_high, require_not_std_vector_t<T_y>* = nullptr, require_std_vector_t<T_high>* = nullptr>
 inline void check_less(const char* function, const char* name,
                                    const T_y& y, const T_high& high) {
   for (size_t i = 0; i < y.size(); ++i) {
     check_less(function, name, y, high[i]);
   }
 }

}  // namespace math
}  // namespace stan
#endif

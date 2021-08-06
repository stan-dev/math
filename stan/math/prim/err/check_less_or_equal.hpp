#ifndef STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_mat.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>

namespace stan {
namespace math {

  template <typename T_y, typename T_high, require_all_stan_scalar_t<T_y, T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
     check_not_nan(function, name, y);
     check_not_nan(function, "higher", high);
     if (!(y <= high)) {
       [&]() STAN_COLD_PATH {
         std::stringstream msg;
         msg << ", but must be less than ";
         msg << high;
         std::string msg_str(msg.str());
         throw_domain_error(function, name, y, "is ", msg_str.c_str());
       }();
     }
  }

  template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
   require_matrix_t<T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
     auto&& high_arr = as_array_or_scalar(to_ref(value_of_rec(high)));
     check_not_nan(function, name, y);
     check_not_nan(function, "higher", high_arr);
     for (Eigen::Index j = 0; j < high_arr.cols(); ++j) {
       for (Eigen::Index i = 0; i < high_arr.rows(); ++i) {
       if (!(y <= high_arr.coeffRef(i, j))) {
         [&high_arr, y, name, function, i, j]() STAN_COLD_PATH {
           std::stringstream msg;
           msg << ", but must be less than ";
           msg << high_arr.coeff(i, j);
           std::string msg_str(msg.str());
           throw_domain_error(function, name, y, "is ", msg_str.c_str());
         }();
       }
     }}
  }

template <typename T_y, typename T_high, require_matrix_t<T_y>* = nullptr,
  require_stan_scalar_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
     auto&& y_arr = as_array_or_scalar(to_ref(value_of_rec(y)));
     check_not_nan(function, name, y_arr);
     check_not_nan(function, "higher", high);
     for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
       for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
       if (!(y_arr.coeffRef(i, j) <= high)) {
         [&y_arr, high, name, function, i, j]() STAN_COLD_PATH {
           std::stringstream msg;
           msg << ", but must be less than ";
           msg << high;
           std::string msg_str(msg.str());
           throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                  msg_str.c_str());
         }();
       }
     }}
  }

  template <typename T_y, typename T_high, require_all_matrix_t<T_y, T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
     check_matching_dims(function, name, y, "higher", high);
     auto&& y_arr = to_ref(as_array_or_scalar(value_of_rec(y)));
     auto&& high_arr = to_ref(as_array_or_scalar(value_of_rec(high)));
     check_not_nan(function, name, y_arr);
     check_not_nan(function, "higher", high_arr);
     for (Eigen::Index j = 0; j < y_arr.cols(); ++j) {
       for (Eigen::Index i = 0; i < y_arr.rows(); ++i) {
       if (!(y_arr.coeffRef(i, j) <= high_arr.coeffRef(i, j))) {
         [&y_arr, &high_arr, name, function, i, j]() STAN_COLD_PATH {
           std::stringstream msg;
           msg << ", but must be less than ";
           msg << high_arr.coeff(i, j);
           std::string msg_str(msg.str());
           throw_domain_error_mat(function, name, y_arr, i, j, "is ",
                                  msg_str.c_str());
         }();
       }
     }}
  }

  template <typename T_y, typename T_high, require_all_std_vector_t<T_y, T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
    check_matching_sizes(function, name, y, "higher", high);
    std::string iter_name{name};
    for (size_t i = 0; i < y.size(); ++i) {
      check_less_or_equal(function, (iter_name + "[" + std::to_string(i + stan::error_index::value) + "]").c_str(), y[i], high[i]);
    }
  }

  template <typename T_y, typename T_high, require_std_vector_t<T_y>* = nullptr,
   require_not_std_vector_t<T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
    std::string iter_name{name};
    for (size_t i = 0; i < y.size(); ++i) {
      check_less_or_equal(function, (iter_name + "[" + std::to_string(i + stan::error_index::value) + "]").c_str(), y[i], high);
    }
  }

  template <typename T_y, typename T_high, require_not_std_vector_t<T_y>* = nullptr, require_std_vector_t<T_high>* = nullptr>
  inline void check_less_or_equal(const char* function, const char* name,
                                     const T_y& y, const T_high& high) {
    for (size_t i = 0; i < high.size(); ++i) {
      check_less_or_equal(function, name, y, high[i]);
    }
  }
}  // namespace math
}  // namespace stan
#endif

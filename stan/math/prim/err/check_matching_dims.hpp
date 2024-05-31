#ifndef STAN_MATH_PRIM_ERR_CHECK_MATCHING_DIMS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_MATCHING_DIMS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dims.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the two containers have the same dimensions.
 * @tparam T1 type of the first container
 * @tparam T2 type of the second container
 * @param function name of function (for error messages)
 * @param name1 variable name for the first container (for error messages)
 * @param y1 first container to test
 * @param name2 variable name for the second container (for error messages)
 * @param y2 second container to test
 * @throw <code>std::invalid_argument</code> if the dimensions of the
 *    containers do not match
 */
template <typename T1, typename T2, require_all_not_matrix_t<T1, T2>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T1, T2>* = nullptr>
inline void check_matching_dims(const char* function, const char* name1,
                                const T1& y1, const char* name2, const T2& y2) {
  std::vector<int> y1_d = dims(y1);
  std::vector<int> y2_d = dims(y2);
  auto error_throw = [&]() STAN_COLD_PATH {
    std::ostringstream y1s;
    if (y1_d.size() > 0) {
      y1s << y1_d[0];
      for (int i = 1; i < y1_d.size(); i++) {
        y1s << ", " << y1_d[i];
      }
    }
    std::ostringstream msg;
    msg << ") and " << name2 << " (";
    if (y2_d.size() > 0) {
      msg << y2_d[0];
      for (int i = 1; i < y2_d.size(); i++) {
        msg << ", " << y2_d[i];
      }
    }
    msg << ") must match in size";
    std::string msg_str(msg.str());
    invalid_argument(function, name1, y1s.str(), "(", msg_str.c_str());
  };
  if (y1_d.size() != y2_d.size()) {
    error_throw();
  } else {
    for (int i = 0; i < y1_d.size(); i++) {
      if (y1_d[i] != y2_d[i]) {
        error_throw();
      }
    }
  }
}

/**
 * Check if two matrices have the same row and column dimensions.
 * @tparam T1 Either an Eigen type or a `var_value` with underlying Eigen type.
 * @tparam T2 Either an Eigen type or a `var_value` with underlying Eigen type.
 * @param function name of function (for error messages)
 * @param name1 variable name for the first container (for error messages)
 * @param y1 first matrix to test
 * @param name2 variable name for the second container (for error messages)
 * @param y2 second matrix to test
 * @throw <code>std::invalid_argument</code> if the dimensions of the
 *    containers do not match
 */
template <
    typename T1, typename T2,
    require_any_t<conjunction<is_matrix<T1>, is_matrix<T2>>,
                  conjunction<is_prim_or_rev_kernel_expression<T1>,
                              is_prim_or_rev_kernel_expression<T2>>>* = nullptr,
    require_any_not_stan_scalar_t<T1, T2>* = nullptr>
inline void check_matching_dims(const char* function, const char* name1,
                                const T1& y1, const char* name2, const T2& y2) {
  if (y1.rows() != y2.rows() || y1.cols() != y2.cols()) {
    [&]() STAN_COLD_PATH {
      std::ostringstream y1_err;
      std::ostringstream msg_str;
      y1_err << "(" << y1.rows() << ", " << y1.cols() << ") and ";
      msg_str << " (" << y2.rows() << ", " << y2.cols()
              << ") must match in size";
      invalid_argument(function, name1, name2,
                       std::string(y1_err.str()).c_str(),
                       std::string(msg_str.str()).c_str());
    }();
  }
}

/**
 * Check if two matrices have the same row and column dimensions.
 * @tparam T1 Either an Eigen type, a `var_value` with underlying Eigen type, or
 * scalar.
 * @tparam T2 Either an Eigen type, a `var_value` with underlying Eigen type, or
 * scalar.
 * @param function name of function (for error messages)
 * @param name1 variable name for the first container (for error messages)
 * @param y1 first argument to test
 * @param name2 variable name for the second container (for error messages)
 * @param y2 second argument to test
 * @throw <code>std::invalid_argument</code> if the dimensions of the
 *    containers do not match
 */
template <typename T1, typename T2, require_any_matrix_t<T1, T2>* = nullptr,
          require_any_stan_scalar_t<T1, T2>* = nullptr>
inline void check_matching_dims(const char* function, const char* name1,
                                const T1& y1, const char* name2, const T2& y2) {
  std::string y1_err("");
  std::string msg_str("Tried Checking the dimensions of a matrix vs a scalar");
  invalid_argument(function, name1, y1_err, "", msg_str.c_str());
}

/**
 * Check if the two matrices are of the same size.
 * This function checks the runtime sizes and can also check the static
 * sizes as well. For example, a 4x1 matrix is not the same as a vector
 * with 4 elements.
 * @tparam check_compile Whether to check the static sizes
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 * @param function name of function (for error messages)
 * @param name1 variable name for the first matrix (for error messages)
 * @param y1 first matrix to test
 * @param name2 variable name for the second matrix (for error messages)
 * @param y2 second matrix to test
 * @throw <code>std::invalid_argument</code> if the dimensions of the matrices
 *    do not match
 */
template <bool check_compile, typename Mat1, typename Mat2,
          typename = require_all_eigen_t<Mat1, Mat2>>
inline void check_matching_dims(const char* function, const char* name1,
                                const Mat1& y1, const char* name2,
                                const Mat2& y2) {
  if (check_compile
      && (static_cast<int>(Mat1::RowsAtCompileTime)
              != static_cast<int>(Mat2::RowsAtCompileTime)
          || static_cast<int>(Mat1::ColsAtCompileTime)
                 != static_cast<int>(Mat2::ColsAtCompileTime))) {
    [&]() STAN_COLD_PATH {
      std::ostringstream msg;
      msg << "Static rows and cols of " << name1 << " and " << name2
          << " must match in size.";
      std::string msg_str(msg.str());
      invalid_argument(function, msg_str.c_str(), "", "");
    }();
  }
  check_matching_dims(function, name1, y1, name2, y2);
}

}  // namespace math
}  // namespace stan
#endif

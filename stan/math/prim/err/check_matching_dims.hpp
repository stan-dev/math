#ifndef STAN_MATH_PRIM_ERR_CHECK_MATCHING_DIMS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_MATCHING_DIMS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the two matrices are of the same size.
 * This function checks the runtime sizes only.
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 * @param function name of function (for error messages)
 * @param name1 variable name for the first matrix (for error messages)
 * @param y1 first matrix to test
 * @param name2 variable name for the second matrix (for error messages)
 * @param y2 second matrix to test
 * @throw <code>std::invalid_argument</code> if the dimensions of the
 *    matrices do not match
 */
template <typename Mat1, typename Mat2,
          typename = require_all_eigen_t<Mat1, Mat2>>
inline void check_matching_dims(const char* function, const char* name1,
                                const Mat1& y1, const char* name2,
                                const Mat2& y2) {
  check_size_match(function, "Rows of ", name1, y1.rows(), "rows of ", name2,
                   y2.rows());
  check_size_match(function, "Columns of ", name1, y1.cols(), "columns of ",
                   name2, y2.cols());
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
    std::ostringstream msg;
    msg << "Static rows and cols of " << name1 << " and " << name2
        << " must match in size.";
    std::string msg_str(msg.str());
    invalid_argument(function, msg_str.c_str(), "", "");
  }
  check_matching_dims(function, name1, y1, name2, y2);
}

}  // namespace math
}  // namespace stan
#endif

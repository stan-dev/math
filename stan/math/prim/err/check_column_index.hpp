#ifndef STAN_MATH_PRIM_ERR_CHECK_COLUMN_INDEX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_COLUMN_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/out_of_range.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the specified index is a valid column of the matrix.
 * By default this is a 1-indexed check (as opposed to
 * 0-indexed). Behavior can be changed by setting
 * <code>stan::error_index::value</code>. This function will
 * throw an <code>std::out_of_range</code> exception if
 * the index is out of bounds.
 * @tparam T_y Type of matrix
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y matrix to test
 * @param i column index to check
 * @throw <code>std::out_of_range</code> if index is an invalid column
 */
template <typename T_y,
          require_any_t<is_matrix<T_y>,
                        is_prim_or_rev_kernel_expression<T_y>>* = nullptr>
inline void check_column_index(const char* function, const char* name,
                               const T_y& y, size_t i) {
  if (!(i >= stan::error_index::value
        && i < static_cast<size_t>(y.cols()) + stan::error_index::value)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << " for columns of " << name;
      std::string msg_str(msg.str());
      out_of_range(function, y.cols(), i, msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif

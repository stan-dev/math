#ifndef STAN_MATH_PRIM_ERR_CHECK_ROW_INDEX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_ROW_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/out_of_range.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the specified index is a valid row of the matrix
 * This check is 1-indexed by default. This behavior can be changed
 * by setting <code>stan::error_index::value</code>.
 * @tparam T Matrix type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y matrix to test
 * @param i row index to check
 * @throw <code>std::out_of_range</code> if the index is out of range.
 */
template <typename T_y, typename = require_matrix_t<T_y>>
inline void check_row_index(const char* function, const char* name,
                            const T_y& y, size_t i) {
  STAN_NO_RANGE_CHECKS_RETURN;
  if (!(i >= stan::error_index::value
        && i < static_cast<size_t>(y.rows()) + stan::error_index::value)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << " for rows of " << name;
      std::string msg_str(msg.str());
      out_of_range(function, y.rows(), i, msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif

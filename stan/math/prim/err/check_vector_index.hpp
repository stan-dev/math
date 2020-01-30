#ifndef STAN_MATH_PRIM_ERR_CHECK_VECTOR_INDEX_HPP
#define STAN_MATH_PRIM_ERR_CHECK_VECTOR_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/out_of_range.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the specified index is a valid element of the row or column vector
 * This check is 1-indexed by default. This behavior can be changed
 * by setting <code>stan::error_index::value</code>.
 * @tparam T Vector type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y vector to test
 * @param i row index to check
 * @throw <code>std::out_of_range</code> if the index is out of range.
 */
template <typename T, typename = require_eigen_vector_t<T>>
inline void check_vector_index(const char* function, const char* name,
                               const T& y, size_t i) {
  if (i >= stan::error_index::value
      && i < static_cast<size_t>(y.size()) + stan::error_index::value) {
    return;
  }

  std::stringstream msg;
  msg << " for size of " << name;
  std::string msg_str(msg.str());
  out_of_range(function, y.rows(), i, msg_str.c_str());
}

}  // namespace math
}  // namespace stan
#endif

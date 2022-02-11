#ifndef STAN_MATH_PRIM_ERR_CHECK_SQUARE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_SQUARE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <sstream>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is square. This check allows 0x0 matrices.
 * @tparam T Type of matrix
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::invalid_argument</code> if the matrix is not square
 */
template <typename T_y,
          require_any_t<is_matrix<T_y>,
                        is_prim_or_rev_kernel_expression<T_y>>* = nullptr>
inline void check_square(const char* function, const char* name, const T_y& y) {
  check_size_match(function, "Expecting a square matrix; rows of ", name,
                   y.rows(), "columns of ", name, y.cols());
}

}  // namespace math
}  // namespace stan
#endif

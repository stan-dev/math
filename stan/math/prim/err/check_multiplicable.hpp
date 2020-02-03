#ifndef STAN_MATH_PRIM_ERR_CHECK_MULTIPLICABLE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_MULTIPLICABLE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_size_match.hpp>

namespace stan {
namespace math {

/**
 * Throw exception if the matrices cannot be multiplied.
 * This requires the number of the columns of the first matrix
 * to match the number of rows of the second matrix.
 *
 * @tparam T1 type of first matrix, vector, or row vector
 * @tparam T2 type of second matrix, vector, or row vector
 * @param function function name
 * @param name1 name of first matrix
 * @param y1 first matrix
 * @param name2 name of second matrix
 * @param y2 second matrix
 * @throw <code>std::invalid_argument</code> if the matrices are not
 *   multiplicable
 */
template <typename T1, typename T2>
inline void check_multiplicable(const char* function, const char* name1,
                                const T1& y1, const char* name2, const T2& y2) {
  check_size_match(function, "Columns of ", name1, y1.cols(), "Rows of ", name2,
                   y2.rows());
}
}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_OPENCL_ERR_CHECK_MATCHING_DIMS_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_MATCHING_DIMS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

namespace stan {
namespace math {
/** \ingroup error_checks_opencl
 * Check if two <code>matrix_cl</code>s have the same dimensions.
 *
 * @param function Function name (for error messages)
 * @param name1 Variable name for the first matrix (for error messages)
 * @param y1 First <code>matrix_cl</code>
 * @param name2 Variable name for the second matrix (for error messages)
 * @param y2 Second <code>matrix_cl</code>
 *
 * @throw <code>std::invalid_argument</code>
 * if the dimensions of the matrices do not match
 */
template <typename T>
inline void check_matching_dims(const char* function, const char* name1,
                                const matrix_cl<T>& y1, const char* name2,
                                const matrix_cl<T>& y2) {
  check_size_match(function, "Rows of ", name1, y1.rows(), "rows of ", name2,
                   y2.rows());
  check_size_match(function, "Columns of ", name1, y1.cols(), "columns of ",
                   name2, y2.cols());
}

}  // namespace math
}  // namespace stan
#endif
#endif

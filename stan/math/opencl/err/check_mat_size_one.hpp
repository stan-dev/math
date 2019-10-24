#ifndef STAN_MATH_OPENCL_ERR_CHECK_MAT_SIZE_ONE_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_MAT_SIZE_ONE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>

namespace stan {
namespace math {

/**
 * Check if the <code>matrix_cl</code> has a single element.
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param A <code>matrix_cl</code> to test
 *
 * @throw <code>std::domain_error</code> if the <code>matrix_cl</code>
 *    size is not 1
 */
template <typename T>
inline void check_mat_size_one(const char* function, const char* name,
                               const matrix_cl<T>& A) {
  if (A.size() != 1) {
    invalid_argument(function, name, "should have exactly one element.", "");
  }
}

/**
 * Check if the <code>matrix_cl</code> has a single element.
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param A <code>matrix_cl</code> to test
 *
 * @throw <code>std::domain_error</code> if the <code>matrix_cl</code>
 *    size is not 1
 */
template <typename T>
inline void check_mat_not_size_one(const char* function, const char* name,
                                   const matrix_cl<T>& A) {
  if (A.size() == 1) {
    invalid_argument(function, name, "should not have exactly one element.",
                     "");
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif

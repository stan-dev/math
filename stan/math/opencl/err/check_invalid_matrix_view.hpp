#ifndef STAN_MATH_OPENCL_ERR_CHECK_INVALID_MATRIX_VIEW_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_INVALID_MATRIX_VIEW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>

namespace stan {
namespace math {

/**
 * Check if the <code>matrix_cl</code> has an invalid view.
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param A <code>matrix_cl</code> to test
 * @param invalid_view the view that is not allowed
 *
 * @throw <code>std::domain_error</code> if the <code>matrix_cl</code>
 *    size is not 1
 */
template <typename T>
inline void check_invalid_matrix_view(const char* function, const char* name,
                                      const matrix_cl<T>& A,
                                      const matrix_cl_view invalid_view) {
  if (A.view() == invalid_view) {
    invalid_argument(function, name, " has an invalid view.", "");
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif

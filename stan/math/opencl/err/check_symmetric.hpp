#ifndef STAN_MATH_OPENCL_ERR_CHECK_SYMMETRIC_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_SYMMETRIC_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/kernels/check_symmetric.hpp>
#include <vector>

namespace stan {
namespace math {
/** \ingroup error_checks_opencl
 * Check if the <code>matrix_cl</code> is symmetric
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y <code>matrix_cl</code> to test
 *
 * @throw <code>std::domain_error</code> if
 *    the matrix is not symmetric.
 */
template <typename T, typename = require_arithmetic_t<T>>
inline void check_symmetric(const char* function, const char* name,
                            const matrix_cl<T>& y) {
  if (y.size() == 0) {
    return;
  }
  check_square(function, name, y);
  try {
    int symmetric_flag = 1;
    matrix_cl<int> symm_flag(1, 1);
    symm_flag = to_matrix_cl(symmetric_flag);
    opencl_kernels::check_symmetric(cl::NDRange(y.rows(), y.cols()), y,
                                    symm_flag, y.rows(), y.cols(),
                                    math::CONSTRAINT_TOLERANCE);
    symmetric_flag = from_matrix_cl<int>(symm_flag);
    if (!symmetric_flag) {
      throw_domain_error(function, name, "is not symmetric", "");
    }
  } catch (const cl::Error& e) {
    check_opencl_error("symmetric_check", e);
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif

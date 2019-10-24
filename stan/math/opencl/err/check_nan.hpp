#ifndef STAN_MATH_OPENCL_ERR_CHECK_NAN_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_NAN_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernels/check_nan.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

#include <vector>

namespace stan {
namespace math {
/**
 * Check if the <code>matrix_cl</code> has NaN values
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y <code>matrix_cl</code> to test
 *
 * @throw <code>std::domain_error</code> if
 *    any element of the matrix is <code>NaN</code>.
 */
template <typename T, typename = require_floating_point_t<T>>
inline void check_nan(const char* function, const char* name,
                      const matrix_cl<T>& y) {
  if (y.size() == 0) {
    return;
  }
  try {
    int nan_flag = 0;
    matrix_cl<int> nan_chk(1, 1);
    nan_chk = to_matrix_cl(nan_flag);
    opencl_kernels::check_nan(cl::NDRange(y.rows(), y.cols()), y, nan_chk,
                              y.rows(), y.cols());
    nan_flag = from_matrix_cl_error_code(nan_chk);
    if (nan_flag) {
      domain_error(function, name, "has NaN values", "");
    }
  } catch (const cl::Error& e) {
    check_opencl_error("nan_check", e);
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif

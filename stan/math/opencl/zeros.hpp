#ifndef STAN_MATH_OPENCL_ZEROS_HPP
#define STAN_MATH_OPENCL_ZEROS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/zeros.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

#include <cl.hpp>

namespace stan {
namespace math {

/**
 * Stores zeros in the matrix on the OpenCL device.
 * Supports writing zeroes to the lower and upper triangular or
 * the whole matrix.
 *
 * @tparam view Specifies if zeros are assigned to
 * the entire matrix, lower triangular or upper triangular. The
 * value must be of type matrix_cl_view
 */
template <typename T>
template <matrix_cl_view matrix_view>
inline void matrix_cl<T, enable_if_arithmetic<T>>::zeros() try {
  if (size() == 0) {
    return;
  }
  this->view_ = both(this->view_, invert(matrix_view));
  cl::CommandQueue cmdQueue = opencl_context.queue();
  opencl_kernels::zeros(cl::NDRange(this->rows(), this->cols()), *this,
                        this->rows(), this->cols(), matrix_view);
} catch (const cl::Error& e) {
  check_opencl_error("zeros", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif

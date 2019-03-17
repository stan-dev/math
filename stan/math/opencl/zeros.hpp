#ifndef STAN_MATH_OPENCL_ZEROS_HPP
#define STAN_MATH_OPENCL_ZEROS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/kernels/zeros.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Stores zeros in the matrix on the GPU.
 * Supports writing zeroes to the lower and upper triangular or
 * the whole matrix.
 *
 * @tparam triangular_view Specifies if zeros are assigned to
 * the entire matrix, lower triangular or upper triangular. The
 * value must be of type TriangularViewGPU
 */
template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
inline void matrix_cl::zeros() {
  if (size() == 0)
    return;
  try {
    cl::Event zero_event
        = opencl_kernels::zeros(cl::NDRange(this->rows(), this->cols()), this,
                                this->rows(), this->cols(), triangular_view);
    this->events(zero_event);
  } catch (const cl::Error& e) {
    check_opencl_error("zeros", e);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif

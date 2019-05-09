#ifndef STAN_MATH_OPENCL_TRIANGULAR_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/kernels/triangular_transpose.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

#include <CL/cl.hpp>

namespace stan {
namespace math {

/**
 * Copies a lower/upper triangular of a matrix to it's upper/lower.
 *
 * @tparam triangular_map Specifies if the copy is
 * lower-to-upper or upper-to-lower triangular. The value
 * must be of type TriangularMap
 *
 * @throw <code>std::invalid_argument</code> if the matrix is not square.
 *
 */
template <TriangularMapCL triangular_map>
inline void matrix_cl::triangular_transpose() try {
  if (size() == 0 || size() == 1) {
    return;
  }
  check_size_match("triangular_transpose ((OpenCL))",
                   "Expecting a square matrix; rows of ", "A", rows(),
                   "columns of ", "A", cols());

  cl::CommandQueue cmdQueue = opencl_context.queue();
  opencl_kernels::triangular_transpose(cl::NDRange(this->rows(), this->cols()),
                                       *this, this->rows(), this->cols(),
                                       triangular_map);
} catch (const cl::Error& e) {
  check_opencl_error("triangular_transpose", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif

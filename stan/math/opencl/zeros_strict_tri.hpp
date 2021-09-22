#ifndef STAN_MATH_OPENCL_ZEROS_HPP
#define STAN_MATH_OPENCL_ZEROS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/kernels/fill_strict_tri.hpp>
#include <stan/math/prim/meta.hpp>

#include <CL/opencl.hpp>

namespace stan {
namespace math {

/** \ingroup matrix_cl_group
 * Stores zeros in the strict's triangular part (excluding the diagonal)
 * of a matrix on the OpenCL device.
 * Supports writing zeroes to the lower and upper triangular.
 * Throws if used with the Entire matrix_cl_view.
 *
 * @tparam view Specifies if zeros are assigned to
 * the lower triangular or upper triangular. The
 * value must be of type matrix_cl_view
 *
 * @throw <code>std::invalid_argument</code> if the
 * matrix_view parameter is Entire.
 */
template <typename T>
template <matrix_cl_view matrix_view>
inline void matrix_cl<T>::zeros_strict_tri() try {
  if (matrix_view == matrix_cl_view::Entire) {
    invalid_argument(
        "zeros_strict_tri", "matrix_view",
        "matrix_cl_view::Entire is not a valid template parameter value", "");
  }
  if (matrix_view == matrix_cl_view::Diagonal) {
    invalid_argument(
        "zeros_strict_tri", "matrix_view",
        "matrix_cl_view::Diagonal is not a valid template parameter value", "");
  }
  if (size() == 0) {
    return;
  }
  this->view_ = both(this->view_, invert(matrix_view));
  cl::CommandQueue cmdQueue = opencl_context.queue();
  opencl_kernels::fill_strict_tri(cl::NDRange(this->rows(), this->cols()),
                                  *this, 0.0, this->rows(), this->cols(),
                                  matrix_view);
} catch (const cl::Error& e) {
  check_opencl_error("zeros_strict_tri", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif

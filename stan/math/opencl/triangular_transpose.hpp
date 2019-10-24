#ifndef STAN_MATH_OPENCL_TRIANGULAR_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/triangular_transpose.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

#include <cl.hpp>

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
template <typename T>
template <TriangularMapCL triangular_map>
inline void matrix_cl<T, require_arithmetic_t<T>>::triangular_transpose() try {
  if (this->size() == 0 || this->size() == 1) {
    return;
  }
  check_size_match("triangular_transpose ((OpenCL))",
                   "Expecting a square matrix; rows of ", "A", this->rows(),
                   "columns of ", "A", this->cols());

  cl::CommandQueue cmdQueue = opencl_context.queue();
  opencl_kernels::triangular_transpose(cl::NDRange(this->rows(), this->cols()),
                                       *this, this->rows(), this->cols(),
                                       triangular_map);
  this->view_
      = (triangular_map == TriangularMapCL::LowerToUpper
         && !contains_nonzero(this->view_, matrix_cl_view::Lower))
                || (triangular_map == TriangularMapCL::UpperToLower
                    && !contains_nonzero(this->view_, matrix_cl_view::Upper))
            ? matrix_cl_view::Diagonal
            : matrix_cl_view::Entire;
} catch (const cl::Error& e) {
  check_opencl_error("triangular_transpose", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif

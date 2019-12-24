#ifndef STAN_MATH_OPENCL_REV_ZEROS_HPP
#define STAN_MATH_OPENCL_REV_ZEROS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <stan/math/opencl/zeros.hpp>

#include <CL/cl2.hpp>

namespace stan {
namespace math {

/**
 * Stores zeros in the matrix on the OpenCL device.
 * Supports writing zeroes to the lower and upper triangular or
 * the whole matrix.
 *
 * @tparam partial_view Specifies if zeros are assigned to
 * the entire matrix, lower triangular or upper triangular. The
 * value must be of type matrix_cl_view
 */
template <typename T>
template <matrix_cl_view matrix_view>
inline void matrix_cl<T, require_var_t<T>>::zeros() try {
  this->val().template zeros<matrix_view>();
  this->adj().template zeros<matrix_view>();
} catch (const cl::Error& e) {
  check_opencl_error("zeros", e);
}

/**
 * Stores zeros in the stricts triangular part (excluding the diagonal)
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
inline void matrix_cl<T, require_var_t<T>>::zeros_strict_tri() try {
  this->val().template zeros_strict_tri<matrix_view>();
  this->adj().template zeros_strict_tri<matrix_view>();
} catch (const cl::Error& e) {
  check_opencl_error("zeros_strict_tri", e);
}

}  // namespace math
}  // namespace stan

#endif
#endif

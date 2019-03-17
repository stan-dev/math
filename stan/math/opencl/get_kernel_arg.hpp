#ifndef STAN_MATH_OPENCL_GET_KERNEL_ARG_HPP
#define STAN_MATH_OPENCL_GET_KERNEL_ARG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace opencl_kernels {

/**
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @tparam For this general template the function will just return back the
 * value passed in.
 * @param t The type that will be returned.
 * @return the input t.
 */
template <typename T>
inline const T& get_kernel_arg(const T& t) {
  return t;
}

/**
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @param m The matrix to extract the buffer from.
 * @return the input t's buffer.
 */
inline const cl::Buffer& get_kernel_arg(const matrix_cl& m) {
  return m.buffer();
}

/**
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @param m The matrix to extract the buffer from.
 * @return the input t's buffer.
 */
inline const cl::Buffer& get_kernel_arg(matrix_cl* const& m) {
  return m->buffer();
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

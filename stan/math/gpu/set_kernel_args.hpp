#ifndef STAN_MATH_GPU_SET_KERNEL_ARGS_HPP
#define STAN_MATH_GPU_SET_KERNEL_ARGS_HPP
#ifdef STAN_OPENCL
#include <CL/cl.hpp>

namespace stan {
namespace math {
inline void _set_kernel_args(cl::Kernel &k, int i) {}  // terminating function

/**
 * Used in <code>set_kernel_args()</code> to add arguments to an OpenCL kernel.
 *
 * @param kernel An OpenCL kernel.
 * @param i the position of the argument to the OpenCL kernel.
 * @param first_arg The first argument to go into the OpenCL kernel.
 * @param extra_args The remaining arguments to go into the OpenCL kernel.
 * @tparam T the type of <code>first_arg</code>.
 * @tparam Args The types of <code>extra_args</code>.
 * @note Comes from:
 * simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
 */
template <typename T, typename... Args>
inline void _set_kernel_args(cl::Kernel &kernel, int i, const T &first_arg,
                             const Args &... extra_args) {
  kernel.setArg(i, first_arg);
  _set_kernel_args(kernel, i + 1, extra_args...);
}

/**
 * Adds arguments to an OpenCL kernel.
 *
 * @param kernel An OpenCL kernel.
 * @param args The arguments to mote to the OpenCL kernel.
 * @tparam Args The types of <code>extra_args</code>.
 * @note Comes from:
 * simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
 */
template <typename... Args>
inline void set_kernel_args(cl::Kernel &kernel, const Args &... args) {
  _set_kernel_args(kernel, 0, args...);
}
}  // namespace math
}  // namespace stan
#endif
#endif

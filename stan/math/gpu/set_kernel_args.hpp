#ifndef STAN_MATH_GPU_SET_KERNEL_ARGS_HPP
#define STAN_MATH_GPU_SET_KERNEL_ARGS_HPP
#ifdef STAN_OPENCL
#include <CL/cl.hpp>

namespace stan {
namespace math {
inline void _set_kernel_args(cl::Kernel &k, int i) {}  // terminating function

// simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
template <typename T, typename... Args>
inline void _set_kernel_args(cl::Kernel &kernel, int i, const T &first_arg,
                             const Args &... extra_args) {
  kernel.setArg(i, first_arg);
  _set_kernel_args(kernel, i + 1, extra_args...);
}

template <typename... Args>
inline void set_kernel_args(cl::Kernel &kernel, const Args &... args) {
  _set_kernel_args(kernel, 0, args...);
}
}  // namespace math
}  // namespace stan
#endif
#endif

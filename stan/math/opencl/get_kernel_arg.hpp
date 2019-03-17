#ifndef STAN_MATH_OPENCL_GET_KERNEL_ARG_HPP
#define STAN_MATH_OPENCL_GET_KERNEL_ARG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace opencl_kernels {

template <typename T>
inline const T& get_kernel_arg(const T& t) {
  return t;
}

inline const cl::Buffer& get_kernel_arg(const matrix_cl& m) {
  return m.buffer();
}

inline const cl::Buffer& get_kernel_arg(matrix_cl* const& m) {
  return m->buffer();
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_SELECT_EVENTS_HPP
#define STAN_MATH_OPENCL_SELECT_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <CL/cl.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

/**
 * Specialization for non-matrix_cl types for getting event stack from
 * matrix_cls.
 * @tparam A non-matrix_cl type.
 * @return An empty vector.
 */
template <typename T>
inline const std::vector<cl::Event> select_events(
    internal::to_const_matrix_cl_v<T>& t) {
  return std::vector<cl::Event>();
}

/**
 * Gets the event stack for in_buffers
 * @param m A matrix_cl holding an event stack.
 * @return The write event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<in_buffer>(
    const matrix_cl& m) {
  return m.write_events();
}

/**
 * Gets the event stack for read and write buffers
 * @param m A matrix_cl holding an event stack.
 * @return The read event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<out_buffer>(
    const matrix_cl& m) {
  return m.read_write_events();
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

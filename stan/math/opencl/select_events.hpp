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
inline const std::vector<cl::Event> select_events(const T& t) {
  return std::vector<cl::Event>();
}

/**
 * Gets the event stack for read and write buffers
 * @tparam event_buffer Whether the matrix is to be treated as read or write.
 * @param m A matrix_cl holding an event stack.
 * @return Depending on the template type will return either the read or
 * read_write event stacks.
 */
template <typename event_buffer>
inline const std::vector<cl::Event> select_events(const matrix_cl& m) {}


/**
 * Gets the event stack for read_buffers
 * @param m A matrix_cl holding an event stack.
 * @return The write event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<read_buffer>(const matrix_cl& m) {
  return m.write_events();
}

/**
 * Gets the event stack for read and write buffers
 * @param m A matrix_cl holding an event stack.
 * @return The read and write event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<write_buffer>(const matrix_cl& m) {
  return m.read_write_events();
}

/**
 * Gets the event stack for read and write buffers
 * @tparam event_buffer Whether the matrix is to be treated as read or write.
 * @param m A matrix_cl holding an event stack.
 * @return Depending on the template type will return either the read or
 * read_write event stacks.
 */
template <typename event_buffer>
inline const std::vector<cl::Event> select_events(matrix_cl* const& m) {}


/**
 * Gets the event stack for read buffers
 * @param m Pointer to matrix_cl holding an event stack.
 * @return The write event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<read_buffer>(matrix_cl* const& m) {
  return m->write_events();
}

/**
 * Gets the event stack for write buffers
 * @param m Pointer to matrix_cl holding an event stack.
 * @return The read and write event stack.
 */
template <>
inline const std::vector<cl::Event> select_events<write_buffer>(matrix_cl* const& m) {
  return m->read_write_events();
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

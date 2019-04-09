#ifndef STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#define STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>

#include <CL/cl.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

/**
 * Generic template for non-matrix_cl types
 * @tparam T Any non-matrix-cl type will be declared here and not assigned an
 * event.
 * @param t A type that is not a matrix_cl.
 * @param new_event An event to add to an event stack
 */
template <typename T>
inline void assign_event(const T& t, const cl::Event& new_event) {}

/**
 * Generic template for matrix_cl types
 * @tparam buffer_type Whether the matrix_cl should add to it's read or write
 * stack.
 * @param m Matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <typename buffer_type>
inline void assign_event(const matrix_cl& m, const cl::Event& new_event) {}

/**
 * Adds event to matrices appropriate event stack.
 * @param m Matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <>
inline void assign_event<read_buffer>(const matrix_cl& m,
                                      const cl::Event& new_event) {
  m.add_read_event(new_event);
}

/**
 * Adds event to matrices appropriate event stack.
 * @param m Matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <>
inline void assign_event<write_buffer>(const matrix_cl& m,
                                       const cl::Event& new_event) {
  m.add_write_event(new_event);
}

/**
 * Generic template for matrix_cl pointer types
 * @tparam buffer_type Whether the matrix_cl should add to it's read or write
 * stack.
 * @param m Matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <typename buffer_type>
inline void assign_event(const matrix_cl* m, const cl::Event& new_event) {}

/**
 * Adds event to matrices appropriate event stack.
 * @param m Pointer to matrix_cl to add to write event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <>
inline void assign_event<read_buffer>(const matrix_cl* m,
                                      const cl::Event& new_event) {
  m->add_read_event(new_event);
}

/**
 * Adds event to matrices appropriate event stack.
 * @param m Pointer to matrix_cl to add to write event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <>
inline void assign_event<write_buffer>(const matrix_cl* m,
                                       const cl::Event& new_event) {
  m->add_write_event(new_event);
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

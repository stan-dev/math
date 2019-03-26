#ifndef STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#define STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

/**
 * Generic template for non-matrix_cl types
 * @tparam Any non-matrix-cl type will be declared here and not assign an event.
 */
template <typename T>
inline void assign_event(const T& t, cl::Event new_event) {}

/**
 * Adds event to matrices appropriate event stack.
 * @tparam buffer_type Whether the matrix_cl should add to it's read or write
 * stack.
 * @param m Matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <typename buffer_type = write_buffer>
inline void assign_event(const matrix_cl& m, cl::Event new_event) {
  m.add_event<buffer_type::event_type>(new_event);
}

/**
 * Adds event to matrices appropriate event stack.
 * @param m Pointer to matrix to add the event to it's event stack.
 * @param new_event The event to add to the matrices event stack.
 */
template <typename buffer_type = write_buffer>
inline void assign_event(const matrix_cl*& m, cl::Event new_event) {
  m->add_event<buffer_type::event_type>(new_event);
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

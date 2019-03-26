#ifndef STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#define STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

template <typename T>
inline void assign_event(const T& t, cl::Event new_event) {}

template <typename P = read_write_buffer>
inline void assign_event(const matrix_cl& m, cl::Event new_event) {
  m.add_event<P::event_type>(new_event);
}

template <typename P = read_write_buffer>
inline void assign_event(const matrix_cl* & m, cl::Event new_event) {
   m->add_event<P::event_type>(new_event);
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

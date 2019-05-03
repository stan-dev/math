#ifndef STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#define STAN_MATH_OPENCL_ASSIGN_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

#include <CL/cl.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

namespace internal {
template <typename T>
void assign_event(const cl::Event &, to_const_matrix_cl_v<T> &) {}

template <>
void assign_event<in_buffer>(const cl::Event &e, const matrix_cl &m) {
  m.add_read_event(e);
}

template <>
void assign_event<out_buffer>(const cl::Event &e, const matrix_cl &m) {
  m.add_write_event(e);
}

template <>
void assign_event<in_out_buffer>(const cl::Event &e, const matrix_cl &m) {
  m.add_read_event(e);
  m.add_write_event(e);
}
}  // namespace internal

template <typename Arg,
          typename std::enable_if<std::is_same<Arg, cl::Event>::value>::type
              * = nullptr>
inline void assign_events(const Arg &) {}

/**
 * Adds the event to any matrices in the arguments in the event vector specified
 * by the buffer directionality.
 * @tparam Arg Arguments given during kernel creation that specify the kernel
 * signature.
 * @tparam Args Arguments given during kernel creation that specify the kernel
 * signature.
 * @param new_event The cl::Event generated involving the arguments.
 * @param m Arguments to the kernel that may be matrices or not. Non-matrices
 * ignored.
 * @param args Arguments to the kernel that may be matrices or not. Non-matrices
 * ignored.
 */
template <typename Arg, typename... Args>
inline void assign_events(const cl::Event &new_event,
                          internal::to_const_matrix_cl_v<Arg> &m,
                          internal::to_const_matrix_cl_v<Args> &... args) {
  internal::assign_event<Arg>(new_event, m);
  assign_events<Args...>(new_event, args...);
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

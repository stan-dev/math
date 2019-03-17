#ifndef STAN_MATH_OPENCL_GET_EVENT_CL_HPP
#define STAN_MATH_OPENCL_GET_EVENT_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

template <typename T>
inline const std::vector<cl::Event> get_event_cl(const T& t) {
  return std::vector<cl::Event>();
}

inline const std::vector<cl::Event>& get_event_cl(const matrix_cl& m) {
  return m.events();
}

inline const std::vector<cl::Event>& get_event_cl(matrix_cl* const& m) {
  return m->events();
}

inline const std::vector<cl::Event>& get_event_cl(
    const std::vector<cl::Event>& m) {
  return m;
}
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_SELECT_EVENTS_HPP
#define STAN_MATH_OPENCL_SELECT_EVENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

template <typename T>
inline const std::vector<cl::Event> select_events(const T& t) {
  return std::vector<cl::Event>();
}

inline const std::vector<cl::Event>& select_events(const matrix_cl& m) {
  return m.events();
}

inline const std::vector<cl::Event>& select_events(
    const std::vector<cl::Event>& m) {
  return m;
}
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

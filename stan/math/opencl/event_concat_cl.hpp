#ifndef STAN_MATH_OPENCL_EVENT_CONCAT_CL_HPP
#define STAN_MATH_OPENCL_EVENT_CONCAT_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/has_event_stack.hpp>
#include <CL/cl.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

/**
 * Ends the recursion to extract the event stack.
 * @param v1 Events on the OpenCL event stack.
 * @return returns input vector.
 */
inline const std::vector<cl::Event>& event_concat_cl(
    const std::vector<cl::Event>& v1) {
  return v1;
}

/**
 * Passes over an argument without an event stack.
 * @param throwaway_val A value that does not have an event stack.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam T checks if either has an event stack and if not then fails
 * compilation.
 * @tparam Args Types for variadic.
 * @return An empty event vector.
 */
template <typename T, enable_if_no_event_stack<T> = 0, typename... Args>
inline const std::vector<cl::Event> event_concat_cl(const T& throwaway_val,
                                                    const Args... args) {
  return event_concat_cl(args...);
}

/**
 * Gets the event stack from a vector of events and other arguments.
 * @param v1 A event stack to roll up.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam Args Types for variadic.
 * @return Vector of OpenCL events
 */
template <typename... Args>
inline const std::vector<cl::Event> event_concat_cl(
    const std::vector<cl::Event>& v1, const Args... args) {
  std::vector<cl::Event> vec_concat = event_concat_cl(args...);
  vec_concat.insert(vec_concat.end(), v1.begin(), v1.end());
  return vec_concat;
}

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

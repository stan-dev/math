#ifndef STAN_MATH_OPENCL_EVENT_CONCAT_CL_HPP
#define STAN_MATH_OPENCL_EVENT_CONCAT_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Template check for if type has event stack
 * @tparam Type to check.
 * @return enum with false value if type does not have an event stack.
 */
template <typename T>
struct has_event_stack {
  enum { value = false };
};

/**
 * Template check for if type has event stack
 * @tparam Type to check.
 * @return enum with false value if type does not have an event stack.
 */
template <>
struct has_event_stack<matrix_cl> {
  enum { value = true };
};

/**
 * Template check for if type has event stack
 * @tparam Type to check.
 * @return enum with true value if type does not have an event stack.
 */
template <>
struct has_event_stack<std::vector<cl::Event>> {
  enum { value = true };
};

/**
 * Helper template that enables function if type has event stack.
 * @tparam Type to check.
 * @return None, will fail compilation and pass to next template ala SFINAE.
 */
template <typename T>
using enable_if_has_event_stack
    = std::enable_if_t<has_event_stack<T>::value, int>;

/**
 * Helper template that enables function if type does not have event stack.
 * @tparam Type to check.
 * @return None, will fail compilation and pass to next template ala SFINAE.
 */
template <typename T>
using enable_if_no_event_stack
    = std::enable_if_t<!has_event_stack<T>::value, int>;

/**
 * Ends the recurstion to extract the event stack.
 * @return An empty event vector.
 */
inline const void event_concat_cl() {}

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
 * Ends the recursion to extract the event stack.
 * @param A OpenCL matrix holding the events.
 * @return The event stack for the matrix.
 */
inline const std::vector<cl::Event>& event_concat_cl(const matrix_cl& A) {
  return A.events();
}

/**
 * Ends the recurstion to extract the event stack.
 * @param throwaway_val A value that does not have an event stack.
 * @tparam T checks if either has an event stack and if not then fails
 * compilation.
 * @return An empty event vector.
 */
template <typename T, enable_if_no_event_stack<T> = 0>
inline const std::vector<cl::Event> event_concat_cl(const T& throwaway_val) {
  const std::vector<cl::Event> vec_concat;
  return vec_concat;
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

/**
 * Get the event stack from a matrix_cl.
 * @param A matrix holding an event stack.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam Args Types for variadic.
 * @return Vector of OpenCL events
 */
template <typename... Args>
inline const std::vector<cl::Event> event_concat_cl(const matrix_cl& A,
                                                    const Args... args) {
  return event_concat_cl(A.events(), args...);
}

/**
 * Recursion when called from within matrix_cl.
 * @param A A pointer to a matrix_cl type, this is called from within a
 * matrix_cl.
 * @param args variadic arcs passed down to the next recursion.
 * @tparam Args Types for variadic.
 * @return Vector of OpenCL events
 */
template <typename... Args>
inline const std::vector<cl::Event> event_concat_cl(const matrix_cl* const& A,
                                                    const Args... args) {
  return event_concat_cl(A->events(), args...);
}
}  // namespace math
}  // namespace stan

#endif
#endif

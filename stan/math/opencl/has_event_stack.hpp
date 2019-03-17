#ifndef STAN_MATH_OPENCL_HAS_EVENT_STACK_HPP
#define STAN_MATH_OPENCL_HAS_EVENT_STACK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <vector>

namespace stan {
namespace math {
namespace opencl_kernels {

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

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif

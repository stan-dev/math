#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_WRAPPER
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_WRAPPER
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
namespace internal {

/**
 * A wrapper for references. This is used to wrap references when putting them
 * in tuples.
 */
template <typename T>
struct wrapper {
  T x;
  explicit wrapper(T&& x) : x(std::forward<T>(x)) {}
};

template <typename T>
wrapper<T> make_wrapper(T&& x) {
  return wrapper<T>(std::forward<T>(x));
}

}  // namespace internal
/** @}*/
}  // namespace math
}  // namespace stan

#endif

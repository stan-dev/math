#ifndef STAN_MATH_OPENCL_REV_TO_ARENA_HPP
#define STAN_MATH_OPENCL_REV_TO_ARENA_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/rev/arena_type.hpp>

namespace stan {
namespace math {

/**
 * Converts given argument into a type that either has any dynamic allocation on
 * AD stack or schedules its destructor to be called when AD stack memory is
 * recovered.
 *
 * This overload is used for `matrix_cl`. It is converted to a subclass that
 * schedules its destructor to be called when the memory is recovered.
 *
 * @tparam T type of scalar
 * @param a argument
 * @return argument
 */
template <typename T>
arena_t<matrix_cl<T>> to_arena(const matrix_cl<T>& a) {
  arena_t<matrix_cl<T>> res(a.buffer(), a.rows(), a.cols(), a.view());
  for (cl::Event e : a.read_events()) {
    res.add_read_event(e);
  }
  for (cl::Event e : a.write_events()) {
    res.add_write_event(e);
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif

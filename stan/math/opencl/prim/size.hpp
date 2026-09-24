#ifndef STAN_MATH_OPENCL_PRIM_SIZE_HPP
#define STAN_MATH_OPENCL_PRIM_SIZE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <cstdint>

namespace stan {
namespace math {

/**
 * Returns the size (number of the elements) of a `matrix_cl` or
 * `var_value<matrix_cl<T>>`.
 * @param m input to determine size of
 * @return number of elements in m
 */
template <typename T,
          require_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
int64_t size(const T& m) {
  return m.rows() * m.cols();
}

/**
 * Returns the size of a device scalar (`opencl::ScalarCl`), which is one.
 * @tparam T type of the device scalar
 * @return 1
 */
template <typename T, require_scalar_cl_t<T>* = nullptr>
inline constexpr int64_t size(const T& /* x */) noexcept {
  return 1;
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_PRIM_SD_HPP
#define STAN_MATH_OPENCL_PRIM_SD_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/variance.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Returns the unbiased sample standard deviation of the
 * coefficients in the specified std vector, column vector, row vector, or
 * matrix.
 *
 * @tparam T type of the container
 *
 * @param a Specified container.
 * @return Sample standard deviation.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline double sd(const T& a) {
  return sqrt(variance(a));
}

}  // namespace math
}  // namespace stan

#endif
#endif

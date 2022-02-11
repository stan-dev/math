#ifndef STAN_MATH_OPENCL_TO_REF_FOR_OPENCL_HPP
#define STAN_MATH_OPENCL_TO_REF_FOR_OPENCL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/ref_type_for_opencl.hpp>

namespace stan {
namespace math {

/**
 * Converts given Eigen expression into one that can be directly copied to an
 * OpenCL device to create `matrix_cl`. If given expression can be directly
 * copied, this is a no-op that should be optimized away.
 * @tparam T argument type
 * @param a argument
 * @return optionally evaluated argument
 */
template <typename T>
inline ref_type_for_opencl_t<T&&> to_ref_for_opencl(T&& a) {
  return std::forward<T>(a);
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_PRIM_SOFTMAX_HPP
#define STAN_MATH_OPENCL_PRIM_SOFTMAX_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector.
 *
 * @tparam ColVec type of elements in the vector
 * @param[in] v Vector to transform.
 * @return Unit simplex result of the softmax transform of the vector.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto softmax(const T& a) {
  check_vector("softmax (OpenCL)","a",a);
  if (v.size() == 0) {
    return v;
  }
  auto a_eval = eval(a);


}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_PRIM_LOG_SOFTMAX_HPP
#define STAN_MATH_OPENCL_PRIM_LOG_SOFTMAX_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/log_sum_exp.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
inline matrix_cl<double> log_softmax(const T& a) {
  check_nonzero_size("log_softmax (OpenCL)", "x", a);
  return make_holder_cl([](const auto& x) { return x - log_sum_exp(x); },
                        to_ref(a));
}

}  // namespace math
}  // namespace stan

#endif
#endif

#ifndef STAN_MATH_OPENCL_VALUE_TYPE_HPP
#define STAN_MATH_OPENCL_VALUE_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_traits
 * Return the value type of an OpenCL matrix.
 */
template <typename T>
struct value_type<T, require_all_kernel_expressions_and_none_scalar_t<T>> {
  using type = typename std::decay_t<T>::Scalar;
};
}  // namespace stan
#endif
#endif

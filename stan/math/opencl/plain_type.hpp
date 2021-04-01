#ifndef STAN_MATH_OPENCL_PLAIN_TYPE_HPP
#define STAN_MATH_OPENCL_PLAIN_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/prim/meta/plain_type.hpp>

namespace stan {

/**
 * Determines plain (non expression) type associated with \c T. For kernel
 * generator expression it is a type the expression can be evaluated into.
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type<T, require_all_kernel_expressions_and_none_scalar_t<T>> {
  using type = math::matrix_cl<typename std::decay_t<T>::Scalar>;
};

}  // namespace stan

#endif
#endif

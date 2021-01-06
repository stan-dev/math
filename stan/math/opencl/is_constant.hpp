#ifndef STAN_MATH_OPENCL_IS_CONSTANT_HPP
#define STAN_MATH_OPENCL_IS_CONSTANT_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_constant.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

namespace stan {

/** \ingroup type_trait
 * Defines a static member named value and sets it to true
 * if the type of the elements in the provided matrix_cl
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 * @tparam type of the elements in the matrix_cl
 */
template <typename T>
struct is_constant<T, require_all_kernel_expressions_and_none_scalar_t<T>>
    : std::true_type {};

}  // namespace stan

#endif
#endif

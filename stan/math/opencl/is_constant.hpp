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
struct is_constant<math::matrix_cl<T>> : is_constant<std::decay_t<T>> {};

}  // namespace stan

#endif
#endif

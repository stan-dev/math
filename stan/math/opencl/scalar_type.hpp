#ifndef STAN_MATH_OPENCL_SCALAR_TYPE_HPP
#define STAN_MATH_OPENCL_SCALAR_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_traits
 * Return the scalar type of an OpenCL matrix.
 */
template <typename T>
struct scalar_type<T, require_matrix_cl_t<T>> {
  using type = typename scalar_type<typename std::decay_t<T>::Scalar>::type;
};
}  // namespace stan
#endif
#endif

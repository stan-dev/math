#ifndef STAN_MATH_OPENCL_SCALAR_TYPE_HPP
#define STAN_MATH_OPENCL_SCALAR_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct scalar_type<T, std::enable_if_t<is_matrix_cl<T>::value>> {
  using type = typename scalar_type<typename T::Scalar>::type;
};
}  // namespace stan
#endif
#endif

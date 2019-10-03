#ifndef STAN_MATH_OPENCL_VALUE_TYPE_HPP
#define STAN_MATH_OPENCL_VALUE_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct value_type<T, require_matrix_cl_t<T>> {
  using type = typename std::decay_t<T>::Scalar;
};
}  // namespace stan
#endif
#endif

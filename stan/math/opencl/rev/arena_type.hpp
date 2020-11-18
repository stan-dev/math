#ifndef STAN_MATH_OPENCL_REV_ARENA_TYPE_HPP
#define STAN_MATH_OPENCL_REV_ARENA_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/meta/arena_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/rev/arena_matrix_cl.hpp>

namespace stan {
namespace internal {

template <typename T>
struct arena_type_impl<stan::math::matrix_cl<T>> {
  using type = stan::math::arena_matrix_cl<T>;
};

template <typename T>
struct arena_type_impl<
    T, require_all_t<is_kernel_expression_and_not_scalar<T>,
                     bool_constant<!is_matrix_cl<T>::value>,
                     bool_constant<!is_arena_matrix_cl<T>::value>>> {
  using type =
      typename arena_type_impl<stan::math::matrix_cl<value_type_t<T>>>::type;
};

}  // namespace internal
}  // namespace stan
#endif
#endif

#ifndef STAN_MATH_OPENCL_REV_ARENA_TYPE_HPP
#define STAN_MATH_OPENCL_REV_ARENA_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/meta/arena_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>

namespace stan {
namespace internal {

template <typename T>
struct arena_type_impl<T, require_matrix_cl_t<T>> {
  struct type_impl : public stan::math::chainable_alloc, public T {
    using Scalar = typename T::Scalar;
    using type = typename T::type;
    template <typename... Args>
    type_impl(Args&&... args)
        : stan::math::chainable_alloc(), T(std::forward<Args>(args)...) {}
  };
  // workaround restriction of having typedef with the same name as the class
  using type = type_impl;
};

template <typename T>
struct arena_type_impl<T,
                       require_all_t<is_kernel_expression_and_not_scalar<T>,
                                     bool_constant<!is_matrix_cl<T>::value>>> {
  using type = typename arena_type_impl<stan::math::matrix_cl<value_type_t<T>>>::type;
};

}  // namespace internal
}  // namespace stan
#endif
#endif

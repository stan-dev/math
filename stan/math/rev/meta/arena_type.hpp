#ifndef STAN_MATH_REV_META_ARENA_TYPE_HPP
#define STAN_MATH_REV_META_ARENA_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/rev/core/arena_allocator.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/var_value_fwd_declare.hpp>

namespace stan {

namespace internal {
template <typename T, typename = void, typename = void>
struct arena_type_impl {};

template <typename T>
struct arena_type_impl<
    T, require_t<std::is_base_of<stan::math::chainable_alloc, T>>> {
  using type = T;
};

template <typename T>
struct arena_type_impl<T, require_all_t<std::is_trivially_destructible<T>,
                                        bool_constant<!is_eigen<T>::value>>> {
  using type = T;
};

template <typename T, typename Alloc>
struct arena_type_impl<std::vector<T, Alloc>> {
  using T_ad = typename arena_type_impl<std::decay_t<T>>::type;
  using type = std::vector<T_ad, math::arena_allocator<T_ad>>;
};

template <typename T>
struct arena_type_impl<
    T, require_eigen_t<T>,
    std::enable_if_t<T::RowsAtCompileTime == Eigen::Dynamic
                     || T::ColsAtCompileTime == Eigen::Dynamic>> {
  using type = math::arena_matrix<plain_type_t<T>>;
};

template <typename T>
struct arena_type_impl<
    T, require_eigen_t<T>,
    std::enable_if_t<T::RowsAtCompileTime != Eigen::Dynamic
                     && T::ColsAtCompileTime != Eigen::Dynamic>> {
  using type = plain_type_t<T>;
};
}  // namespace internal

/**
 * Determines a type that can be used in place of `T` that does any dynamic
 * allocations on the AD stack. This way resulting types are trivially
 * destructible and can be used in vari classes.
 */
template <typename T>
using arena_t = typename internal::arena_type_impl<std::decay_t<T>>::type;

}  // namespace stan

#endif

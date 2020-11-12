#ifndef STAN_MATH_PRIM_META_IS_DENSE_DYNAMIC_HPP
#define STAN_MATH_PRIM_META_IS_DENSE_DYNAMIC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_eigen_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

namespace internal {
template <typename T, typename = void>
struct is_dense_dynamic_impl : std::false_type {};

template <typename T>
struct is_dense_dynamic_impl<T,
                             require_t<is_eigen_dense_dynamic<std::decay_t<T>>>>
    : std::true_type {};

template <typename T>
struct is_dense_dynamic_impl<T, require_t<is_var<T>>>
    : bool_constant<is_eigen_dense_dynamic<value_type_t<T>>::value> {};
}  // namespace internal

/**
 * Checks whether type T is derived from Eigen::DenseBase and has dynamic rows
 * and columns or is a `var_value<>` whose inner type satisfies the conditions
 * above. If true this will have a static member function named value with
 * a type of true, else value is false.
 * @tparam T Type to check
 * @ingroup type_trait
 */
template <typename T>
using is_dense_dynamic = internal::is_dense_dynamic_impl<std::decay_t<T>>;

STAN_ADD_REQUIRE_UNARY(dense_dynamic, is_dense_dynamic, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(dense_dynamic, is_dense_dynamic,
                           require_eigens_types);

}  // namespace stan

#endif

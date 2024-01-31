#ifndef STAN_MATH_PRIM_META_IS_EIGEN_INDEX_VIEW_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_INDEX_VIEW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <type_traits>

namespace stan {
namespace internal {
template <typename T>
struct is_eigen_index_view : std::false_type {};

template <typename ArgType, typename RowIndices, typename ColIndices>
struct is_eigen_index_view<Eigen::IndexedView<ArgType, RowIndices, ColIndices>>
    : std::true_type {};

}  // namespace internal

/**
 * Checks whether type T is derived from Eigen::DenseBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `DenseBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_index_view : internal::is_eigen_index_view<std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(eigen_index_view, is_eigen_index_view,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_index_view, is_eigen_index_view,
                           require_eigens_types);


} // namespace stan

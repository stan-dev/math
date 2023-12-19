#ifndef STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>
#include <tuple>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Promote a scalar to another scalar type
 *
 * @tparam PromotionScalar scalar type of output.
 * @tparam UnPromotedType input type. `UnPromotedType` must be constructible
 * from `PromotionScalar`
 * @param x input scalar to be promoted to `PromotionScalar` type
 */
template <typename PromotionScalar, typename UnPromotedType,
          require_constructible_t<PromotionScalar, UnPromotedType>* = nullptr,
          require_not_same_t<PromotionScalar, UnPromotedType>* = nullptr,
          require_all_not_tuple_t<PromotionScalar, UnPromotedType>* = nullptr>
inline constexpr auto promote_scalar(UnPromotedType&& x) {
  return PromotionScalar(std::forward<UnPromotedType>(x));
}

/**
 * No-op overload when promoting a type's scalar to the type it already has.
 *
 * @tparam PromotionScalar scalar type of output.
 * @tparam UnPromotedType input type. `UnPromotedType`'s `scalar_type` must be
 * equal to `PromotionScalar`
 * @param x input
 */
template <
    typename PromotionScalar, typename UnPromotedType,
    require_same_t<PromotionScalar, scalar_type_t<UnPromotedType>>* = nullptr>
inline constexpr auto promote_scalar(UnPromotedType&& x) noexcept {
  return std::forward<UnPromotedType>(x);
}

/**
 * Promote the scalar type of an eigen matrix to the requested type.
 *
 * @tparam PromotionScalar scalar type of output.
 * @tparam UnPromotedType input type. The `PromotionScalar` type must be
 * constructible from `UnPromotedType`'s `scalar_type`
 * @param x input
 */
template <typename PromotionScalar, typename UnPromotedType,
          require_eigen_t<UnPromotedType>* = nullptr,
          require_not_same_t<PromotionScalar,
                             value_type_t<UnPromotedType>>* = nullptr>
inline auto promote_scalar(UnPromotedType&& x) {
  return x.template cast<PromotionScalar>();
}

// Forward decl for iterating over tuples used in std::vector<tuple>
template <typename PromotionScalars, typename UnPromotedTypes,
          require_all_tuple_t<PromotionScalars, UnPromotedTypes>* = nullptr,
          require_not_same_t<PromotionScalars, UnPromotedTypes>* = nullptr>
inline constexpr promote_scalar_t<PromotionScalars, UnPromotedTypes>
promote_scalar(UnPromotedTypes&& x);

/**
 * Promote the scalar type of an standard vector to the requested type.
 *
 * @tparam PromotionScalar scalar type of output.
 * @tparam UnPromotedType input type. The `PromotionScalar` type must be
 * constructible from `UnPromotedType`'s `scalar_type`
 * @param x input
 */
template <typename PromotionScalar, typename UnPromotedType,
          require_std_vector_t<UnPromotedType>* = nullptr,
          require_not_same_t<PromotionScalar,
                             scalar_type_t<UnPromotedType>>* = nullptr>
inline auto promote_scalar(UnPromotedType&& x) {
  const auto x_size = x.size();
  promote_scalar_t<PromotionScalar, UnPromotedType> ret(x_size);
  for (size_t i = 0; i < x_size; ++i) {
    ret[i] = promote_scalar<PromotionScalar>(x[i]);
  }
  return ret;
}

/**
 * Promote the scalar type of a tuples elements to the requested types.
 *
 * @tparam PromotionScalars A tuple of scalar types that is the same size as the
 * tuple of `UnPromotedTypes`.
 * @tparam UnPromotedTypes tuple input. Each `PromotionScalars` element must be
 * constructible from it's associated element of `UnPromotedTypes` `scalar_type`
 * @param x input
 */
template <typename PromotionScalars, typename UnPromotedTypes,
          require_all_tuple_t<PromotionScalars, UnPromotedTypes>*,
          require_not_same_t<PromotionScalars, UnPromotedTypes>*>
inline constexpr promote_scalar_t<PromotionScalars, UnPromotedTypes>
promote_scalar(UnPromotedTypes&& x) {
  return index_apply<std::tuple_size<std::decay_t<UnPromotedTypes>>::value>(
      [&x](auto... Is) {
        return std::make_tuple(
            promote_scalar<std::decay_t<decltype(std::get<Is>(
                std::declval<PromotionScalars>()))>>(std::get<Is>(x))...);
      });
}

}  // namespace math
}  // namespace stan

#endif

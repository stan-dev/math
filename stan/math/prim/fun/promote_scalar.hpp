#ifndef STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>
#include <tuple>
#include <type_traits>

namespace stan {
namespace math {


template <typename PromotionScalars, typename UnPromotedTypes>
inline constexpr auto promote_scalar(UnPromotedTypes&& x) {
  if constexpr (std::is_same_v<PromotionScalars, scalar_type_t<UnPromotedTypes>>) {
    return std::forward<UnPromotedTypes>(x);
  } else if constexpr (is_tuple_v<PromotionScalars> && is_tuple_v<UnPromotedTypes>) {
    return index_apply<std::tuple_size<std::decay_t<UnPromotedTypes>>::value>(
      [&x](auto... Is) {
        return std::make_tuple(
            promote_scalar<std::decay_t<decltype(std::get<Is>(
                std::declval<PromotionScalars>()))>>(std::get<Is>(x))...);
      });
  } else if constexpr (is_tuple_v<UnPromotedTypes>) {
    return stan::math::apply([](auto&&... args) {
      return std::make_tuple(promote_scalar<PromotionScalars>(std::forward<decltype(args)>(args))...);
    }, std::forward<UnPromotedTypes>(x));
  } else if constexpr (is_std_vector_v<UnPromotedTypes>) {
    const auto x_size = x.size();
    promote_scalar_t<PromotionScalars, UnPromotedTypes> ret(x_size);
    for (size_t i = 0; i < x_size; ++i) {
      ret[i] = promote_scalar<PromotionScalars>(x[i]);
    }
    return ret;
  } else if constexpr (is_eigen_v<UnPromotedTypes>) {
    return std::forward<UnPromotedTypes>(x).template cast<PromotionScalars>();
  } else if constexpr (is_stan_scalar_v<UnPromotedTypes>) {
    return PromotionScalars(std::forward<UnPromotedTypes>(x));
  } else {
    static_assert(1, "Missed type in promote_scalar!");
  }
}


}  // namespace math
}  // namespace stan

#endif

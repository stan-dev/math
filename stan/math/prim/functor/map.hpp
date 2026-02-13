#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan::math {
template <typename F, typename T>
inline decltype(auto) map(F&& f, T&& arg) {
  if constexpr (is_tuple_v<T>) {
    return stan::math::apply(
        [](auto&& f, auto&&... args) {
          return make_holder_tuple(
              std::forward<decltype(f)>(f)(
              std::forward<decltype(args)>(args))...);
        },
        std::forward<T>(arg), std::forward<F>(f));
  } else {
    return std::forward<F>(f)(std::forward<T>(arg));
  }
}
}

#endif  // STAN_MATH_PRIM_FUNCTOR_MAP_HPP

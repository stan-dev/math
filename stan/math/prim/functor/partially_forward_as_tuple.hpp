#ifndef STAN_MATH_PRIM_FUNCTOR_PARTIALLY_FORWARD_AS_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_PARTIALLY_FORWARD_AS_TUPLE_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>


namespace stan {
namespace math {

template <typename T>
struct deduce_cvr {
  using type =
      std::conditional_t<std::is_rvalue_reference_v<T>,
      std::decay_t<T>, T&&>;
};

template <typename T>
using deduce_cvr_t = typename deduce_cvr<T>::type;

template<typename... Types>
constexpr std::tuple<deduce_cvr_t<Types&&>...> partially_forward_as_tuple( Types&&... args ) noexcept {
    return {std::forward<Types>(args)...};
}
}
}
#endif
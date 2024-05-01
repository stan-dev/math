#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_AT_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_AT_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <boost/mp11.hpp>
#include <tuple>

namespace stan {
namespace math {
namespace internal {

template <typename F, typename T, typename... TArgs>
struct require_same_result_type_impl {
  using type = require_all_same_t<
    decltype(std::declval<F>()(std::declval<T>())),
    decltype(std::declval<F>()(std::declval<TArgs>()))...
  >;
};

template <typename F, typename T, typename... TArgs>
struct require_same_result_type_impl<F, std::tuple<T, TArgs...>> {
  using type = require_all_same_t<
    decltype(std::declval<F>()(std::declval<T>())),
    decltype(std::declval<F>()(std::declval<TArgs>()))...
  >;
};

template <typename F, typename T>
using require_same_result_type_t =
      typename require_same_result_type_impl<F, T>::type;
}

/**
 * Call a functor f at a runtime-specified index of a tuple. This requires that
 * the return type of the functor is identical for every tuple element.
 *
 * @tparam F Type of functor
 * @tparam TupleT Type of tuple containing arguments
 * @param func Functor callable
 * @param t Tuple of arguments
 * @param element Element of tuple to apply functor to
 */
  template <typename F, typename TupleT,
            internal::require_same_result_type_t<F, TupleT>* = nullptr>
  decltype(auto) apply_at(F&& func, size_t element, TupleT&& t) {
    constexpr size_t tuple_size = std::tuple_size<std::decay_t<TupleT>>{};
    return boost::mp11::mp_with_index<tuple_size>(
      element, [&](auto I) -> decltype(auto) {
        return func(std::forward<decltype(std::get<I>(t))>(std::get<I>(t)));
      });
  }
}  // namespace math
}  // namespace stan

#endif

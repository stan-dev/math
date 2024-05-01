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
using require_same_result_type_t = typename require_same_result_type_impl<F, T>::type;
}

/*
 * Call the functor f with the tuple of arguments t, like:
 *
 * f(std::get<0>(t), std::get<1>(t), ...)
 *
 * TODO: replace this with implementation in C++ std when C++17 is available
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @tparam PreArgs Parameter pack of arguments before the tuple
 * @param f functor callable
 * @param t tuple of arguments
 * @param pre_args parameter pack of arguments to place before elements in
 * tuple.
 */
  template <typename F, typename TupleT,
            internal::require_same_result_type_t<F, TupleT>* = nullptr>
  decltype(auto) apply_at(F&& func, size_t element, TupleT&& v) {
    constexpr size_t tuple_size = std::tuple_size<std::decay_t<TupleT>>{};
    return boost::mp11::mp_with_index<tuple_size>(
      element, [&](auto I) -> decltype(auto) {
        return func(std::forward<decltype(std::get<I>(v))>(std::get<I>(v)));
      });
  }
}  // namespace math
}  // namespace stan

#endif

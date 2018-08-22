#ifndef STAN_MATH_PRIM_SCAL_FUN_CALL_ALL_ARGUMENT_COMBOS_HPP
#define STAN_MATH_PRIM_SCAL_FUN_CALL_ALL_ARGUMENT_COMBOS_HPP

#include <tuple>

namespace stan {
namespace math {

template <typename F>
auto call_all_argument_combos(F f) {
  return std::make_tuple(f());
}

template <typename F, typename... Ts_first_arg, std::size_t... I,
          typename... T_tail>
auto call_all_argument_combos_impl(
    F f, const std::tuple<Ts_first_arg...>& first_arg_tuple,
    std::index_sequence<I...>, const T_tail&... tail);

template <typename F, typename... Ts_first_arg, typename... T_tail>
auto call_all_argument_combos(
    F f, const std::tuple<Ts_first_arg...>& first_arg_tuple,
    const T_tail&... tail) {
  return call_all_argument_combos_impl(
      f, first_arg_tuple, std::make_index_sequence<sizeof...(Ts_first_arg)>{},
      tail...);
}

template <typename F, typename... Ts_first_arg, std::size_t... I,
          typename... T_tail>
auto call_all_argument_combos_impl(
    F f, const std::tuple<Ts_first_arg...>& first_arg_tuple,
    std::index_sequence<I...>, const T_tail&... tail) {
  return std::tuple_cat(call_all_argument_combos(
      [&first_arg_tuple, &f](const auto&... inner_args) {
        return f(std::get<I>(first_arg_tuple), inner_args...);
      },
      tail...)...);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_CALL_ALL_ARGUMENT_COMBOS_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_CALL_ALL_ARGUMENT_COMBOS_HPP

#include <tuple>

namespace stan {
namespace math {

template <typename F>
auto call_all_argument_combos(F f) {
  return std::make_tuple(f());
}

/**
 * There needs to be a forward declaration here because call_all_argument_combos
 * and call_all_argument_combos_impl call each other recursively
 */
template <typename F, typename... Ts_first_arg, std::size_t... I,
          typename... T_tail>
auto call_all_argument_combos_impl(
    F f, const std::tuple<Ts_first_arg...>& first_arg_tuple,
    std::index_sequence<I...>, const T_tail&... tail);

/**
 * call_all_argument_combos is a utility for calling a function with all
 * combinations of input arguments specified
 *
 * The arguments to call_all_argument_combos are:
 * 1. a functor f with M input arguments (M >= 1) and and a non-void return type
 * 2. M tuples (t1, t2, ...), one for each argument f
 *
 * f is called each combinations of the elements of the tuples, and returns a
 * flat tuple with all the return values arranged in row-major order (last index
 * moves first).
 *
 * If f takes two arguments, this would look like:
 *
 * f(std::get<0>(t1), std::get<0>(t2))
 * f(std::get<0>(t1), std::get<1>(t2))
 * ...
 * f(std::get<1>(t1), std::get<0>(t2))
 * ...
 *
 * @tparam F type of functor
 * @tparam ...Ts_first_arg types of tuple of first argument
 * @tparam ...T_tail Tuple types of the rest of the arguments
 * @param f functor
 * @param first_arg_tuple Tuple of values for first argument
 * @param T_tail Tuples of values for the trailing arguments
 */
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

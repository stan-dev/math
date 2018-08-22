#ifndef STAN_MATH_PRIM_ARR_FUN_PROMOTE_DOUBLE_TO_T_HPP
#define STAN_MATH_PRIM_ARR_FUN_PROMOTE_DOUBLE_TO_T_HPP

#include <tuple>
#include <vector>
#include <iostream>

namespace stan {
namespace math {

template <typename T, typename... Pargs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output) {
  return output;
}

template <typename T, typename R, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output,
                                   const std::vector<R>& arg,
                                   const Targs&... args);
template <typename T, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output,
                                   const double& arg, const Targs&... args);
template <typename T, typename R, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output, const R& arg,
                                   const Targs&... args);

template <typename T, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output,
                                   const std::vector<double>& arg,
                                   const Targs&... args) {
  std::vector<T> output_arg;
  output_arg.reserve(arg.size());
  for (int i = 0; i < arg.size(); ++i)
    output_arg.push_back(arg[i]);
  return promote_double_to_T_impl_impl<T>(
      std::tuple_cat(output, std::make_tuple(output_arg)), args...);
}

template <typename T, typename R, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output,
                                   const std::vector<R>& arg,
                                   const Targs&... args) {
  return promote_double_to_T_impl_impl<T>(output, args...);
}

template <typename T, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output,
                                   const double& arg, const Targs&... args) {
  return promote_double_to_T_impl_impl<T>(
      std::tuple_cat(output, std::make_tuple(T(arg))), args...);
}

template <typename T, typename R, typename... Pargs, typename... Targs>
auto promote_double_to_T_impl_impl(std::tuple<Pargs...> output, const R& arg,
                                   const Targs&... args) {
  return promote_double_to_T_impl_impl<T>(output, args...);
}

template <typename T, std::size_t... I, typename... Targs>
auto promote_double_to_T_impl(std::index_sequence<I...>,
                              const std::tuple<Targs...>& args) {
  return promote_double_to_T_impl_impl<T>(std::tuple(), std::get<I>(args)...);
}

template <typename T, typename... Targs>
auto promote_double_to_T(const std::tuple<Targs...>& args) {
  return promote_double_to_T_impl<T>(
      std::make_index_sequence<sizeof...(Targs)>{}, args);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_FUN_A_HPP
#define STAN_MATH_PRIM_FUN_A_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <iostream>

namespace stan {
namespace math {
namespace internal {
  template <template <class> class Check>
  constexpr size_t type_count(size_t x) {
    return x;
  }
  template <template <class> class Check, typename T, typename... TArgs>
  constexpr size_t type_count(size_t x) {
    return Check<T>::value ? type_count<Check, TArgs...>(x + 1) : x;
  }
  template <template <class> class Check, typename... TArgs>
  constexpr size_t type_count() {
      return type_count<Check, TArgs...>(0);
}

  template <std::size_t O, std::size_t ... Is>
  constexpr std::index_sequence<(O + Is)...>
    add_offset(std::index_sequence<Is...>) {
      return {};
  }

  template<typename Tuple, std::size_t... Ints>
  auto subset_tuple(Tuple&& tuple, std::index_sequence<Ints...>) {
  return std::forward_as_tuple(std::get<Ints>(std::forward<Tuple>(tuple))...);
  }

  template <class F, typename Tuple, size_t... Is>
  auto iter_apply_impl(Tuple t, F f, std::index_sequence<Is...>) {
      return std::make_tuple(
          f(std::get<Is>(t))...
      );
  }

  template <class F, typename... Args>
  auto iter_apply(const std::tuple<Args...>& t, F f) {
      return iter_apply_impl(
          t, f, std::make_index_sequence<sizeof...(Args)>{});
  }

  template <typename T>
  using is_stan_type = disjunction<is_stan_scalar<T>, is_container<T>>;
  template <typename T>
  using is_arith = is_constant_all<T>;
}

template <typename ScalarT, typename ValTupleT, typename ValFun, typename GradFunT,
          require_arithmetic_t<ScalarT>* = nullptr>
auto user_gradients_impl(const ValTupleT& val_tuple, const ValFun& val_fun,
                         const GradFunT& grad_fun_tuple) {
  auto rtn = math::apply([&](auto&&... args) { return val_fun(args...); }, val_tuple);
  decltype(auto) grad_apply = [&](auto&& f) { return math::apply([&](auto&&... args) { return f(args...); }, val_tuple); };
  auto grd = internal::iter_apply(grad_fun_tuple, grad_apply);
  std::cout << rtn << "\n" << std::get<0>(grd) << std::endl;
}


template <typename... TArgs>
auto user_gradients(const TArgs&... args) {

  // Get the number of input variables/parameters
  constexpr size_t pos = internal::type_count<internal::is_stan_type,
                                              TArgs...>();
  constexpr size_t arg_size = sizeof...(TArgs);

  // Create index sequences for all inputs
  constexpr auto input_indices = std::make_index_sequence<pos>{};
  constexpr size_t value_fun_indices = pos;
  constexpr auto grad_fun_indices = internal::add_offset<pos+1>(
                          std::make_index_sequence<arg_size-pos-1>{});

  // Pass all inputs to a tuple
  decltype(auto) args_tuple = std::forward_as_tuple(args...);
  using TupleT = decltype(args_tuple);


  // Check whether all inputs are primitive
  constexpr size_t arith_pos = internal::type_count<internal::is_arith, TArgs...>();
  constexpr bool constant_only = pos == arith_pos;

  decltype(auto) val_tuple =
    internal::subset_tuple(std::forward<TupleT>(args_tuple), input_indices);
  decltype(auto) val_fun = std::get<value_fun_indices>(args_tuple);
  decltype(auto) grad_fun_tuple =
    internal::subset_tuple(std::forward<TupleT>(args_tuple), grad_fun_indices);

  std::cout << constant_only << std::endl;

  return user_gradients_impl<double>(val_tuple, val_fun, grad_fun_tuple);
}

}  // namespace math
}  // namespace stan

#endif

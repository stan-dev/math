#ifndef STAN_MATH_PRIM_FUN_A_HPP
#define STAN_MATH_PRIM_FUN_A_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <iostream>

namespace stan {
namespace math {
namespace internal {
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


template <typename ArgsTupleT, typename ValFunT, typename GradFunTupleT>
auto user_gradients(const ArgsTupleT& args_tuple, const ValFunT& val_fun,
                    const GradFunTupleT& gradfun_tuple) {
  using rtn_t = return_type_t<ArgsTupleT>;
  return user_gradients_impl<scalar_type_t<rtn_t>>(args_tuple, val_fun, gradfun_tuple);
}

}  // namespace math
}  // namespace stan

#endif

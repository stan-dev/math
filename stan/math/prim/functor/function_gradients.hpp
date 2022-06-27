#ifndef STAN_MATH_PRIM_FUNCTOR_FUNCTION_GRADIENTS_HPP
#define STAN_MATH_PRIM_FUNCTOR_FUNCTION_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename ReturnT, typename ArgsTupleT, typename ValFunT,
          typename RevGradFunT, typename FwdGradFunT,
          require_st_arithmetic<ReturnT>* = nullptr>
decltype(auto) function_gradients_impl(
    ArgsTupleT&& args_tuple, ValFunT&& val_fun,
    RevGradFunT&& rev_grad_fun_tuple, FwdGradFunT&& fwd_grad_fun_tuple) {
  return make_holder(
      [](auto&& fun, auto&& tuple_arg) {
        return to_ref(
            math::apply([&](auto&&... args) { return fun(args...); },
                        std::forward<decltype(tuple_arg)>(tuple_arg)));
      },
      std::forward<ValFunT>(val_fun), std::forward<ArgsTupleT>(args_tuple));
}

template <typename ArgsTupleT, typename ValFunT, typename RevGradFunT,
          typename FwdGradFunT>
decltype(auto) function_gradients(ArgsTupleT&& args_tuple,
                                          ValFunT&& val_fun,
                                          RevGradFunT&& rev_grad_fun_tuple,
                                          FwdGradFunT&& fwd_grad_fun_tuple) {
  using scalar_rtn_t = scalar_type_t<return_type_t<ArgsTupleT>>;
  return function_gradients_impl<scalar_rtn_t>(
      std::forward<ArgsTupleT>(args_tuple), std::forward<ValFunT>(val_fun),
      std::forward<RevGradFunT>(rev_grad_fun_tuple),
      std::forward<FwdGradFunT>(fwd_grad_fun_tuple));
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {

template <typename ScalarT, typename ArgsTupleT,
          typename ValFun, typename GradFunT,
          require_var_t<ScalarT>* = nullptr>
auto user_gradients_impl(ArgsTupleT&& args_tuple, ValFun&& val_fun,
                         GradFunT&& grad_fun_tuple) {
  auto prim_tuple = map_tuple([&](auto&& arg) {
    return to_arena(value_of(arg));
  }, std::forward<ArgsTupleT>(args_tuple));
  auto arena_tuple = map_tuple([&](auto&& arg) {
    return to_arena(arg);
  }, std::forward<ArgsTupleT>(args_tuple));
  auto rtn = math::apply([&](auto&&... args) {
    return val_fun(args...);
  }, std::forward<decltype(prim_tuple)>(prim_tuple));
  
  return make_callback_var(rtn,
    [grad_fun_tuple, arena_tuple, prim_tuple](auto& vi) mutable {
    walk_tuples([&](auto&& f, auto&& arg) {
      using arg_t = plain_type_t<decltype(arg)>;
      if (!is_constant_all<arg_t>::value) {
        forward_as<promote_scalar_t<var, arg_t>>(arg).adj()
          += vi.adj()
              * math::apply([&](auto&&... args) {
                  return f(args...);
                }, std::forward<decltype(prim_tuple)>(prim_tuple));
      }
    }, std::forward<GradFunT>(grad_fun_tuple),
       std::forward<decltype(arena_tuple)>(arena_tuple));
  });
}

}  // namespace math
}  // namespace stan
#endif

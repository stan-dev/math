#ifndef STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {

template <typename ScalarT, typename ValTupleT,
          typename ValFun, typename GradFunT,
          require_var_t<ScalarT>* = nullptr>
auto user_gradients_impl(const ValTupleT& val_tuple, const ValFun& val_fun,
                         const GradFunT& grad_fun_tuple) {
  auto prim_tuple = map_tuple([&](auto&& arg) {
      return to_arena(value_of(arg)); }, val_tuple);
  auto arena_tuple = map_tuple([&](auto&& arg) {
      return to_arena(arg); }, val_tuple);
  auto rtn = math::apply([&](auto&&... args) {
      return val_fun(args...); }, prim_tuple);
  
  return make_callback_var(
      rtn, [grad_fun_tuple, arena_tuple, prim_tuple](auto& vi) mutable {
    walk_tuples([&](auto&& f, auto&& arg) {
      if (!is_constant_all<decltype(arg)>::value) {
        forward_as<promote_scalar_t<var, decltype(arg)>>(arg).adj()
          += vi.adj() * math::apply([&](auto&&... args) {
                                        return f(args...);
                                        }, prim_tuple);
      }
    }, grad_fun_tuple, arena_tuple);
  });
}

}  // namespace math
}  // namespace stan
#endif

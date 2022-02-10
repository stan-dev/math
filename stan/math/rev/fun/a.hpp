#ifndef STAN_MATH_REV_FUN_A_HPP
#define STAN_MATH_REV_FUN_A_HPP

#include <stan/math/prim/fun/a.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {

template <typename ScalarT, typename ValTupleT, typename ValFun, typename GradFunT,
          require_var_t<ScalarT>* = nullptr>
auto user_gradients_impl(const ValTupleT& val_tuple, const ValFun& val_fun,
                         const GradFunT& grad_fun_tuple) {
  auto prim_tuple = internal::iter_apply([&](auto&& arg) { return value_of(arg); }, val_tuple);
  auto arena_tuple = internal::iter_apply([&](auto&& arg) { return to_arena(arg); }, val_tuple);
  auto rtn = math::apply([&](auto&&... args) { return val_fun(args...); }, prim_tuple);
  
  return make_callback_var(rtn, [grad_fun_tuple, arena_tuple](auto& vi) mutable {
    for(int i = 0; i < std::tuple_size<ValTupleT>::value; i++) {
      using elem_t = decltype(std::get<i>(arena_tuple));
      if (!is_constant_all<elem_t>::value) {
        forward_as<promote_scalar_t<var, elem_t>>(std::get<i>(arena_tuple)).adj()
          += vi.adj() * math::apply([&](auto&&... args) { return std::get<i>(grad_fun_tuple)(args...); }, arena_tuple);
      }
    }
  });
}

}  // namespace math
}  // namespace stan
#endif

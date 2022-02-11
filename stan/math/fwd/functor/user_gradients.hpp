#ifndef STAN_MATH_FWD_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_FWD_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/fwd/core/fvar.hpp>

namespace stan {
namespace math {

template <typename ScalarT, typename ArgsTupleT,
          typename ValFun, typename GradFunT,
          require_fvar_t<ScalarT>* = nullptr>
auto user_gradients_impl(ArgsTupleT&& args_tuple,
                         ValFun&& val_fun,
                         GradFunT&& grad_fun_tuple) {
  auto val_tuple = map_tuple([&](auto&& arg) {
    return value_of(arg);
  }, std::forward<ArgsTupleT>(args_tuple));
  auto rtn = math::apply([&](auto&&... args) {
    return val_fun(args...);
  }, std::forward<decltype(val_tuple)>(val_tuple));
  
  decltype(rtn) d_(0);
  
  walk_tuples([&](auto&& f, auto&& arg) {
    using arg_t = plain_type_t<decltype(arg)>;
    if (!is_constant_all<arg_t>::value) {
      d_ += forward_as<promote_scalar_t<ScalarT, arg_t>>(arg).d()
            * math::apply([&](auto&&... args) { return f(args...); },
                          std::forward<decltype(val_tuple)>(val_tuple));
    }
  }, std::forward<GradFunT>(grad_fun_tuple),
     std::forward<ArgsTupleT>(args_tuple));

  return ScalarT(rtn, d_);
}

}  // namespace math
}  // namespace stan
#endif

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
auto user_gradients_impl(const ArgsTupleT& args_tuple, const ValFun& val_fun,
                         const GradFunT& grad_fun_tuple) {
  auto val_tuple = map_tuple([&](auto&& arg) {
      return value_of(arg); }, args_tuple);
  auto rtn = math::apply([&](auto&&... args) {
      return val_fun(args...); }, val_tuple);
  
  decltype(rtn) d_(0);
  
  walk_tuples([&](auto&& f, auto&& arg) {
    if (!is_constant_all<decltype(arg)>::value) {
      d_ +=
        forward_as<promote_scalar_t<ScalarT, decltype(arg)>>(arg).d()
          * math::apply([&](auto&&... args) {
                                      return f(args...);
                                      }, val_tuple);
    }
  }, grad_fun_tuple, args_tuple);

  return ScalarT(rtn, d_);
}

}  // namespace math
}  // namespace stan
#endif

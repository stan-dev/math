#ifndef STAN_MATH_PRIM_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_PRIM_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename ScalarT, typename ValTupleT,
          typename ValFun, typename GradFunT,
          require_arithmetic_t<ScalarT>* = nullptr>
auto user_gradients_impl(ValTupleT&& val_tuple, ValFun&& val_fun,
                         GradFunT&& grad_fun_tuple) {
  return math::apply([&](auto&&... args) {
    return val_fun(args...); }, std::forward<ValTupleT>(val_tuple));
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

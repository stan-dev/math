#ifndef STAN_MATH_FWD_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_FWD_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename T, require_stan_scalar_t<T>* = nullptr>
constexpr T initialize_grad(T&& rtn_val) {
  return 0;
}

template <typename T, require_eigen_t<T>* = nullptr>
plain_type_t<T> initialize_grad(T&& rtn_eigen) {
  return plain_type_t<T>::Zero(rtn_eigen.rows(), rtn_eigen.cols());
}
}  // namespace internal

template <typename ScalarT, typename ArgsTupleT, typename ValFun,
          typename GradFunT, require_st_fvar<ScalarT>* = nullptr>
decltype(auto) user_gradients_impl(ArgsTupleT&& args_tuple, ValFun&& val_fun,
                                   GradFunT&& grad_fun_tuple) {
  decltype(auto) val_tuple
      = map_tuple([&](auto&& arg) { return value_of(arg); },
                  std::forward<ArgsTupleT>(args_tuple));

  decltype(auto) rtn
      = math::apply([&](auto&&... args) { return val_fun(args...); },
                    std::forward<decltype(val_tuple)>(val_tuple));

  auto d_ = internal::initialize_grad(std::forward<decltype(rtn)>(rtn));

  walk_tuples(
      [&](auto&& f, auto&& arg) {
        if (!is_constant_all<decltype(arg)>::value) {
          d_ += math::apply(
              [&](auto&&... args) {
                return f(rtn,
                         forward_as<promote_scalar_t<ScalarT, decltype(arg)>>(arg).d(),
                         args...);
              },
              val_tuple);
        }
      },
      std::forward<GradFunT>(grad_fun_tuple),
      std::forward<ArgsTupleT>(args_tuple));

  return to_fvar(rtn, d_);
}

}  // namespace math
}  // namespace stan
#endif

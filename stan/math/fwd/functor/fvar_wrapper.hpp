#ifndef STAN_MATH_FWD_FUNCTOR_FVAR_WRAPPER_HPP
#define STAN_MATH_FWD_FUNCTOR_FVAR_WRAPPER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/finite_diff_gradient_auto.hpp>
#include <stan/math/fwd/functor/gradient.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <test/unit/math/serializer.hpp>

namespace stan {
namespace math {

template <typename F, typename... TArgs, require_all_not_st_fvar<TArgs...>* = nullptr>
auto fvar_wrapper(const F& func, const TArgs&... args) {
  return func(args...);
}

template <typename F, typename... TArgs, require_any_st_fvar<TArgs...>* = nullptr>
auto fvar_wrapper(const F& func, const TArgs&... args) {
  // Create a tuple of double-only arguments
  auto prim_args = std::make_tuple(stan::math::value_of(args)...);

  using FvarT = return_type_t<TArgs...>;
  using FvarInnerT = typename FvarT::Scalar;

  // Flatten the input arguments to a single column vector
  auto serial_args = stan::math::apply(
    [&](auto&&... tuple_args) {
      return math::to_vector(stan::test::serialize<FvarInnerT>(tuple_args...));
    }, prim_args);

  // Create a 'wrapper' functor which will take the flattened column-vector
  // and transform it to individual arguments which are passed to the
  // user-provided functor
  auto serial_functor = [&](const auto& v) {
    auto ds = stan::test::to_deserializer(v);
    return stan::math::apply([&](auto&&... tuple_args) {
      return func(ds.read(tuple_args)...); }, prim_args);
  };

  FvarInnerT rtn_value;
  Eigen::Matrix<FvarInnerT, -1, 1> grad;
  finite_diff_gradient_auto(serial_functor, serial_args, rtn_value, grad);

  auto ds_grad = stan::test::to_deserializer(grad);
  auto grads = std::make_tuple(ds_grad.read(args)...);
  auto args_tuple = std::forward_as_tuple(args...);
  FvarInnerT rtn_grad = 0;

  math::for_each([&](auto&& fun_grad, auto&& arg){
    if (is_fvar<scalar_type_t<decltype(arg)>>::value) {
      rtn_grad += sum(elt_multiply(fun_grad, forward_as<promote_scalar_t<FvarT, decltype(arg)>>(arg).d()));
    }
    return;
  }, std::forward<decltype(grads)>(grads), std::forward<decltype(args_tuple)>(args_tuple));

  return fvar<FvarInnerT>(rtn_value, rtn_grad);
}
}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_FWD_FUNCTOR_FVAR_WRAPPER_HPP
#define STAN_MATH_FWD_FUNCTOR_FVAR_WRAPPER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/finite_diff_gradient_auto.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <test/unit/math/serializer.hpp>

template <typename T, require_arithmetic_t<T>* = nullptr>
constexpr double grad(const T& x) {
  return 0;
}

template <typename T, require_fvar_t<T>* = nullptr>
auto grad(const T& x) {
  return x.d();
}

template <typename F, typename... TArgs>
auto fvar_wrapper(const F& func, const TArgs&... args) {
  // Create a tuple of double-only arguments
  auto prim_args = std::forward_as_tuple(stan::math::value_of(args)...);

  // Flatten the input arguments to a single column vector
  auto serial_args = stan::math::apply(
    [&](auto&&... tuple_args) {
      return stan::test::serialize_args(tuple_args...);
    }, prim_args);

  // Create a 'wrapper' functor which will take the flattened column-vector
  // and transform it to individual arguments which are passed to the
  // user-provided functor
  auto serial_functor = [&](const auto& v) {
    auto ds = stan::test::to_deserializer(v);
    return stan::math::apply([&](auto&&... tuple_args) {
      return func(ds.read(tuple_args)...); }, prim_args);
  };

  double rtn_value;
  Eigen::VectorXd grad;
  finite_diff_gradient_auto(serial_functor, serial_args, rtn_value, grad);

  auto ds_grad = stan::test::to_deserializer(grad);
  double rtn_grad;

  (void)std::initializer_list<int>{ (rtn_grad += sum(elt_multiply(ds.read(args), grad(args))), 0)... };

  return fvar<double>(rtn_value, rtn_grad);
}

#endif

#ifndef STAN_MATH_FWD_FUNCTOR_FUNCTION_GRADIENTS_HPP
#define STAN_MATH_FWD_FUNCTOR_FUNCTION_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/aggregate_partial.hpp>
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

template <typename ScalarT, typename ArgsTupleT, typename ValFunT,
          typename GradFunT, require_st_fvar<ScalarT>* = nullptr>
decltype(auto) function_gradients_impl(ArgsTupleT&& args_tuple,
                                       ValFunT&& val_fun,
                                       GradFunT&& grad_fun_tuple) {
  decltype(auto) val_tuple
      = map_tuple([&](auto&& arg) { return value_of(arg); },
                  std::forward<ArgsTupleT>(args_tuple));

  decltype(auto) rtn
      = math::apply([&](auto&&... args) { return val_fun(args...); },
                    std::forward<decltype(val_tuple)>(val_tuple));
  using rtn_t = decltype(rtn);

  // g++-4.9 has a bug with using decltype within a lambda for a captured
  // value. Create a tuple of unitialised values of the same type as the
  // function return that will be passed to the grad functor
  plain_type_t<rtn_t> dummy_val;
  auto dummy_tuple = map_tuple([&](auto&& arg) { return dummy_val; },
                               std::forward<ArgsTupleT>(args_tuple));

  auto d_ = internal::initialize_grad(std::forward<rtn_t>(rtn));

  walk_tuples(
      [&](auto&& f, auto&& arg, auto&& dummy) {
        using arg_t = decltype(arg);
        if (!is_constant_all<arg_t>::value) {
          decltype(auto) grad = math::apply(
              [&](auto&&... args) { return f(rtn, args...); }, val_tuple);
          as_array_or_scalar(d_) += aggregate_partial<decltype(dummy)>(
              forward_as<promote_scalar_t<ScalarT, arg_t>>(arg),
              std::forward<decltype(grad)>(grad));
        }
      },
      std::forward<GradFunT>(grad_fun_tuple),
      std::forward<ArgsTupleT>(args_tuple),
      std::forward<decltype(dummy_tuple)>(dummy_tuple));

  return to_fvar(rtn, d_);
}

}  // namespace math
}  // namespace stan
#endif

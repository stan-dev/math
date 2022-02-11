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
  // Extract values from input arguments to use in value
  // and gradient calculations
  auto prim_tuple = map_tuple([&](auto&& arg) {
    return to_arena(value_of(arg));
  }, std::forward<ArgsTupleT>(args_tuple));

  // Copy arguments to arena storage for use in callback
  auto arena_tuple = map_tuple([&](auto&& arg) {
    return to_arena(arg);
  }, std::forward<ArgsTupleT>(args_tuple));

  // Use input values to calculate return value
  auto rtn = math::apply([&](auto&&... args) {
    return val_fun(args...);
  }, std::forward<decltype(prim_tuple)>(prim_tuple));
  
  return make_callback_var(rtn,
    [grad_fun_tuple, arena_tuple, prim_tuple](auto& vi) mutable {
    // Iterate over input arguments, applying the respective gradient function
    // with the tuple of extracted primitive values
    walk_tuples([&](auto&& f, auto&& arg) {
      using arg_t = plain_type_t<decltype(arg)>;
      // Only calculate gradients if the input argument is not primitive
      if (!is_constant_all<arg_t>::value) {
        // Need to wrap the argument in a forward_as<var>() so that it will
        // compile with both primitive and var inputs
        forward_as<promote_scalar_t<var, arg_t>>(arg).adj()
          += vi.adj()
            // Use the relevant gradient function with the tuple of primitive
            // arguments
            * math::apply([&](auto&&... args) {
                return f(args...);
              }, std::forward<decltype(prim_tuple)>(prim_tuple));
      }
    },
    std::forward<GradFunT>(grad_fun_tuple),
    std::forward<decltype(arena_tuple)>(arena_tuple)
    );
  });
}

}  // namespace math
}  // namespace stan
#endif

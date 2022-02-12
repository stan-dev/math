#ifndef STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {
namespace internal {
template <typename T,
          require_arena_matrix_t<T>* = nullptr>
decltype(auto) arena_val(T&& arg) {
  return arg.val();
}
template <typename T,
          require_not_arena_matrix_t<T>* = nullptr>
decltype(auto) arena_val(T&& arg) {
  return arg;
}
}

template <typename ReturnT, typename ArgsTupleT,
          typename ValFun, typename GradFunT,
          require_st_var<ReturnT>* = nullptr>
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

  decltype(auto) f = [&](auto&&... args) {
    return val_fun(internal::arena_val(std::forward<decltype(args)>(args))...);
  };

  // Have to declare the math::apply functor separately to avoid compiler
  // errors using lambda expressions in decltype call
  using val_t = decltype(
    math::apply(f, std::forward<decltype(prim_tuple)>(prim_tuple))
  );
  using ret_type = std::conditional_t<is_any_var_matrix<ArgsTupleT>::value,
                                      return_var_matrix_t<val_t, ArgsTupleT>,
                                      promote_scalar_t<var, val_t>>;
  // Use input values to calculate return value
  arena_t<ret_type> rtn = math::apply(f, std::forward<decltype(prim_tuple)>(prim_tuple));
  
  reverse_pass_callback([grad_fun_tuple, arena_tuple, prim_tuple, rtn]() mutable {
    // Iterate over input arguments, applying the respective gradient function
    // with the tuple of extracted primitive values
    walk_tuples([&](auto&& f, auto&& arg) {
      using arg_t = plain_type_t<decltype(arg)>;
      // Only calculate gradients if the input argument is not primitive
      if (!is_constant_all<arg_t>::value) {
        // Need to wrap the argument in a forward_as<var>() so that it will
        // compile with both primitive and var inputs
        forward_as<promote_scalar_t<var, arg_t>>(arg).adj() += 
          // Use the relevant gradient function with the tuple of primitive
          // arguments
          math::apply([&](auto&&... args) {
              return f(
                rtn.adj(),
                internal::arena_val(std::forward<decltype(args)>(args))...
              );
            }, std::forward<decltype(prim_tuple)>(prim_tuple));
      }
    },
    std::forward<GradFunT>(grad_fun_tuple),
    std::forward<decltype(arena_tuple)>(arena_tuple)
    );
  });
  return ret_type(rtn);
}

}  // namespace math
}  // namespace stan
#endif

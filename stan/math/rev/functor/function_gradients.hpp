#ifndef STAN_MATH_REV_FUNCTOR_FUNCTION_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_FUNCTION_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Helper function for agnostically using arena<T> and scalar
 * types as inputs.
 *
 * Specialisation for non-arena types, reduces to a no-op
 *
 * @tparam T Type of non-arena matrix input
 * @param arg Non-arena matrix input
 * @return Unmodified input
 */
template <typename T, require_not_arena_matrix_t<T>* = nullptr>
inline T arena_val(T&& arg) {
  return std::forward<T>(arg);
}

/**
 * Helper function for agnostically using arena<T> and scalar
 * types as inputs.
 *
 * Specialisation for arena types
 *
 * @tparam T Type of arena matrix input
 * @param arg Arena matrix input
 * @return Value within arena matrix
 */
template <typename T, require_arena_matrix_t<T>* = nullptr>
inline decltype(auto) arena_val(T&& arg) {
  return arg.val();
}
}  // namespace internal

template <bool FwdGradients, typename ReturnT, typename ArgsTupleT,
          typename ValFunT, typename RevGradFunT, typename FwdGradFunT,
          require_st_var<ReturnT>* = nullptr>
decltype(auto) function_gradients_impl(ArgsTupleT&& args_tuple, ValFunT&& val_fun,
                                     RevGradFunT&& rev_grad_fun_tuple,
                                     FwdGradFunT&& fwd_grad_fun_tuple) {
  // Extract values from input arguments to use in value
  // and gradient calculations
  decltype(auto) prim_tuple
      = map_tuple([&](auto&& arg) { return to_arena(value_of(arg)); },
                  std::forward<ArgsTupleT>(args_tuple));

  // Copy arguments to arena storage for use in callback
  decltype(auto) arena_tuple
      = map_tuple([&](auto&& arg) { return to_arena(arg); },
                  std::forward<ArgsTupleT>(args_tuple));

  // Have to wrap the holder functor to avoid compiler
  // errors when using lambda expressions in decltype call
  decltype(auto) holder_wrapper = [&](auto&& fun, auto&& tuple_arg) {
    return make_holder(
        [](auto&& h_f, auto&& h_t) {
          return math::apply(
              [&](auto&&... args) {
                return h_f(
                    internal::arena_val(std::forward<decltype(args)>(args))...);
              },
              std::forward<decltype(h_t)>(h_t));
        },
        std::forward<decltype(fun)>(fun),
        std::forward<decltype(tuple_arg)>(tuple_arg));
  };

  // Get primitive return type of function, used for assessing the need for a
  // var<Matrix> return type
  using val_t = decltype(holder_wrapper(std::forward<ValFunT>(val_fun),
                           std::forward<decltype(prim_tuple)>(prim_tuple)));

  // Assess whether the return type should be var<double>, Matrix<var>,
  // or var<Matrix>
  using ret_type = std::conditional_t<is_any_var_matrix<ArgsTupleT>::value,
                                      return_var_matrix_t<val_t, ArgsTupleT>,
                                      promote_scalar_t<var, val_t>>;

  // Use input values to calculate return value
  arena_t<ret_type> rtn = holder_wrapper(std::forward<ValFunT>(val_fun),
                            std::forward<decltype(prim_tuple)>(prim_tuple));
  if (unlikely(math::size(rtn) == 0)) {
    return ret_type(rtn);
  }

  reverse_pass_callback(
      [rev_grad_fun_tuple, arena_tuple, prim_tuple, rtn]() mutable {
        // Iterate over input arguments, applying the respective gradient
        // function with the tuple of extracted primitive values
        walk_tuples(
            [&](auto&& grad_funs, auto&& arg) {
              // Only calculate gradients if the input argument is not primitive
              if (!is_constant_all<decltype(arg)>::value) {
                // Need to wrap the argument in a forward_as<var>() so that it
                // will compile with both primitive and var inputs
                forward_as<promote_scalar_t<var, decltype(arg)>>(arg).adj() +=
                    // Use the relevant gradient function with the tuple of
                    // primitive arguments
                    math::apply(
                        [&](auto&&... args) {
                          return grad_funs(rtn.val_op(), rtn.adj_op(),
                                   internal::arena_val(
                                       std::forward<decltype(args)>(args))...);
                        },
                        prim_tuple);
              }
            },
            std::forward<RevGradFunT>(rev_grad_fun_tuple),
            std::forward<decltype(arena_tuple)>(arena_tuple));
      });
  return ret_type(rtn);
}

}  // namespace math
}  // namespace stan
#endif

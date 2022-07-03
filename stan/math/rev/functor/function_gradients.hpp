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

/**
 * @brief
 *
 * @tparam FwdGradients
 * @tparam ScalarReturnT
 * @tparam ArgsTupleT
 * @tparam ValFunT
 * @tparam RevGradFunT
 * @tparam FwdGradFunT
 * @param input_args_tuple
 * @param value_functor
 * @param rev_grad_functors_tuple
 * @param fwd_grad_functors_tuple
 * @return decltype(auto)
 */
template <bool FwdGradients, typename ScalarReturnT, typename ArgsTupleT,
          typename ValFunT, typename SharedArgsFunT, typename RevGradFunT,
          typename FwdGradFunT, require_st_var<ScalarReturnT>* = nullptr>
decltype(auto) function_gradients_impl(ArgsTupleT&& input_args_tuple,
                                       ValFunT&& value_functor,
                                       SharedArgsFunT&& shared_args_functor,
                                       RevGradFunT&& rev_grad_functors_tuple,
                                       FwdGradFunT&& fwd_grad_functors_tuple) {
  // Extract values from input arguments to use in value
  // and gradient calculations
  decltype(auto) prim_values_tuple
      = map_tuple([&](auto&& input_arg) { return to_arena(value_of(input_arg)); },
                  std::forward<ArgsTupleT>(input_args_tuple));

  using prim_tuple_t = decltype(prim_values_tuple);

  // Copy arguments to arena storage for use in callback
  decltype(auto) arena_args_tuple
      = map_tuple([&](auto&& input_arg) { return to_arena(input_arg); },
                  std::forward<ArgsTupleT>(input_args_tuple));

  // Have to wrap the holder functor to avoid compiler
  // errors when using lambda expressions in decltype call
  decltype(auto) holder_wrapper = [&](auto&& functor, auto&& tuple_args) {
    return make_holder(
        [](auto&& wrapped_functor, auto&& wrapped_tuple) {
          return math::apply(
              [&](auto&&... args) {
                return wrapped_functor(
                    internal::arena_val(std::forward<decltype(args)>(args))...);
              },
              std::forward<decltype(wrapped_tuple)>(wrapped_tuple));
        },
        std::forward<decltype(functor)>(functor),
        std::forward<decltype(tuple_args)>(tuple_args));
  };

  // Get primitive return type of function, used for assessing the need for a
  // var<Matrix> return type
  using val_t = decltype(holder_wrapper(std::forward<ValFunT>(value_functor),
                           std::forward<prim_tuple_t>(prim_values_tuple)));

  // Assess whether the return type should be var<double>, Matrix<var>,
  // or var<Matrix>
  using ret_type = std::conditional_t<is_any_var_matrix<ArgsTupleT>::value,
                                      return_var_matrix_t<val_t, ArgsTupleT>,
                                      promote_scalar_t<var, val_t>>;

  // Use input values to calculate return value
  arena_t<ret_type> rtn_value
    = holder_wrapper(std::forward<ValFunT>(value_functor),
                     std::forward<prim_tuple_t>(prim_values_tuple));

  if (unlikely(math::size(rtn_value) == 0)) {
    return ret_type(rtn_value);
  }

  reverse_pass_callback(
      [rev_grad_functors_tuple, arena_args_tuple, prim_values_tuple,
        shared_args_functor, rtn_value]() mutable {
        decltype(auto) shared_args_tuple = math::apply(
          [&](auto&&... prim_args) {
            return shared_args_functor(
              rtn_value.val_op(), rtn_value.adj_op(),
              internal::arena_val(
                  std::forward<decltype(prim_args)>(prim_args))...)
          }, prim_values_tuple
        );
        // Iterate over input arguments, applying the respective gradient
        // function with the tuple of extracted primitive values
        walk_tuples(
            [&](auto&& grad_funs, auto&& input_arg) {
              // Only calculate gradients if the input argument is not primitive
              if (!is_constant_all<decltype(input_arg)>::value) {
                // Need to wrap the argument in a forward_as<var>() so that it
                // will compile with both primitive and var inputs
                forward_as<promote_scalar_t<var, decltype(input_arg)>>(input_arg).adj() +=
                    // Use the relevant gradient function with the tuple of
                    // primitive arguments
                    math::apply(
                        [&](auto&&... prim_args) {
                          return grad_funs(rtn_value.val_op(), rtn_value.adj_op(),
                                   shared_args_tuple,
                                   internal::arena_val(
                                       std::forward<decltype(prim_args)>(prim_args))...);
                        },
                        prim_values_tuple);
              }
            },
            std::forward<RevGradFunT>(rev_grad_functors_tuple),
            std::forward<decltype(arena_args_tuple)>(arena_args_tuple));
      });
  return ret_type(rtn_value);
}

}  // namespace math
}  // namespace stan
#endif

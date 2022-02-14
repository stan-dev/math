#ifndef STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_USER_GRADIENTS_HPP

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
 * Implementation function for applying a user-defined function and gradients.
 * This specialisation is for var return types, and is compatible with all
 * var, matvar, and matvar input/output types.
 *
 * The provided gradient functions are only evaluated if the respective input
 * is not a primitive type
 *
 * @tparam Scalar type of return value, used for SFINAE
 * @tparam ArgsTupleT Type of arguments tuple to use with function
 * @tparam ValFun Type of user-defined function
 * @tparam GradFunT Type of gradient functions tuple
 * @param args_tuple Tuple of arguments to be evaluated with function
 * @param val_fun Functor to evaluate provided arguments
 * @param grad_fun_tuple Tuple of functors for calculating gradients wrt to
 * each input.
 * @return Result of applying functor to arguments within provided tuple
 */
template <typename ReturnT, typename ArgsTupleT, typename ValFun,
          typename GradFunT, require_st_var<ReturnT>* = nullptr>
auto user_gradients_impl(ArgsTupleT&& args_tuple, ValFun&& val_fun,
                         GradFunT&& grad_fun_tuple) {
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
  decltype(auto) f = [&](auto&& fun, auto&& tuple_arg) {
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
  using val_t = decltype(f(std::forward<ValFun>(val_fun),
                           std::forward<decltype(prim_tuple)>(prim_tuple)));

  // Assess whether the return type should be var<double>, Matrix<var>,
  // or var<Matrix>
  using ret_type = std::conditional_t<is_any_var_matrix<ArgsTupleT>::value,
                                      return_var_matrix_t<val_t, ArgsTupleT>,
                                      promote_scalar_t<var, val_t>>;

  // Use input values to calculate return value
  arena_t<ret_type> rtn = f(std::forward<ValFun>(val_fun),
                            std::forward<decltype(prim_tuple)>(prim_tuple));

  reverse_pass_callback(
      [grad_fun_tuple, arena_tuple, prim_tuple, rtn]() mutable {
        // Iterate over input arguments, applying the respective gradient
        // function with the tuple of extracted primitive values
        walk_tuples(
            [&](auto&& f, auto&& arg) {
              // Only calculate gradients if the input argument is not primitive
              if (!is_constant_all<decltype(arg)>::value) {
                decltype(auto) grad = math::apply([&](auto&&... args) {
                          return f(rtn.val(),
                                   internal::arena_val(
                                       std::forward<decltype(args)>(args))...);
                        },
                        prim_tuple);
                // Need to wrap the argument in a forward_as<var>() so that it
                // will compile with both primitive and var inputs
                as_array_or_scalar(forward_as<promote_scalar_t<var, decltype(arg)>>(arg)).adj() +=
                    // Use the relevant gradient function with the tuple of
                    // primitive arguments
                      aggregate_partial<plain_type_t<decltype(arg)>>(
                        rtn.adj(),
                        std::forward<decltype(grad)>(grad)
                        );
              }
            },
            std::forward<GradFunT>(grad_fun_tuple),
            std::forward<decltype(arena_tuple)>(arena_tuple));
      });
  return ret_type(rtn);
}

}  // namespace math
}  // namespace stan
#endif

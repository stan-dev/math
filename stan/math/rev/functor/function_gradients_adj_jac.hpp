#ifndef STAN_MATH_REV_FUNCTOR_FUNCTION_GRADIENTS_ADJ_JAC_HPP
#define STAN_MATH_REV_FUNCTOR_FUNCTION_GRADIENTS_ADJ_JAC_HPP

#include <stan/math/rev/fun/aggregate_partial.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/function_gradients.hpp>

namespace stan {
namespace math {


template <typename ReturnT, typename ArgsTupleT, typename ValFunT,
          typename RevGradFunT, typename FwdGradFunT,
          require_st_var<ReturnT>* = nullptr>
auto function_gradients_adj_jac_impl(ArgsTupleT&& args_tuple,
                             ValFunT&& val_fun,
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
  using val_t = decltype(f(std::forward<ValFunT>(val_fun),
                           std::forward<decltype(prim_tuple)>(prim_tuple)));

  // Assess whether the return type should be var<double>, Matrix<var>,
  // or var<Matrix>
  using ret_type = std::conditional_t<is_any_var_matrix<ArgsTupleT>::value,
                                      return_var_matrix_t<val_t, ArgsTupleT>,
                                      promote_scalar_t<var, val_t>>;

  // Use input values to calculate return value
  arena_t<ret_type> rtn = f(std::forward<ValFunT>(val_fun),
                            std::forward<decltype(prim_tuple)>(prim_tuple));

  reverse_pass_callback(
      [rev_grad_fun_tuple, arena_tuple, prim_tuple, rtn]() mutable {
        // Iterate over input arguments, applying the respective gradient
        // function with the tuple of extracted primitive values
        walk_tuples(
            [&](auto&& f, auto&& arg) {
              // Only calculate gradients if the input argument is not primitive
              if (!is_constant_all<decltype(arg)>::value) {
                // Need to wrap the argument in a forward_as<var>() so that it
                // will compile with both primitive and var inputs
                forward_as<promote_scalar_t<var, decltype(arg)>>(arg)
                    .adj()
                    +=
                    // Use the relevant gradient function with the tuple of
                    // primitive arguments
                    math::apply(
                    [&](auto&&... args) {
                      return f(rtn.val().eval(),
                               rtn.adj().eval(),
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

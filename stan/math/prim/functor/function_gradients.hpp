#ifndef STAN_MATH_PRIM_FUNCTOR_FUNCTION_GRADIENTS_HPP
#define STAN_MATH_PRIM_FUNCTOR_FUNCTION_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
/**
 * Implementation functor for applying user-specified gradients for a function
 * with arithmetic inputs, simply returning the function result
 *
 * @tparam FwdGradients Boolean indicating whether user provided functors for
 *          forward-mode
 * @tparam ScalarReturnT Scalar return type of the function
 * @tparam ArgsTupleT Type of tuple containing all input arguments
 * @tparam ValFunT Type of functor for calculating return value
 * @tparam RevGradFunT Type of tuple containing functors for reverse-mode
 *          gradients
 * @tparam FwdGradFunT Type of tuple containing functors for forward-mode
 *          gradients
 * @param input_input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param rev_grad_functors_tuple Tuple containing functors for reverse-mode
 *          gradients
 * @param fwd_grad_functors_tuple Tuple containing functors for forward-mode
 *          gradients
 * @return double or container of double with the values calculated by
 *          the value_functor functor
 */
template <bool FwdGradients, typename ScalarReturnT, typename ArgsTupleT,
          typename ValFunT, typename SharedArgsFunT,
          typename RevGradFunT, typename FwdGradFunT,
          require_st_arithmetic<ScalarReturnT>* = nullptr>
decltype(auto) function_gradients_impl(
    ArgsTupleT&& input_args_tuple, ValFunT&& value_functor,
    SharedArgsFunT&& shared_args_functor,
    RevGradFunT&& rev_grad_functors_tuple,
    FwdGradFunT&& fwd_grad_functors_tuple) {
  // The function return needs to be wrapped in a Holder<T> to allow for
  // returning Eigen expressions
  return make_holder(
      [](auto&& functor, auto&& tuple_args) {
        return math::apply([&](auto&&... input_args) { return functor(input_args...); },
                        std::forward<decltype(tuple_args)>(tuple_args));
      },
      std::forward<ValFunT>(value_functor),
      std::forward<ArgsTupleT>(input_args_tuple));
}

/**
 * Entry point for the function_gradients framework. This function checks the
 * scalar return type of the function and then delegates to the appropriate
 * function_gradients_impl overload
 *
 * @tparam ArgsTupleT Type of tuple containing all input arguments
 * @tparam ValFunT Type of functor for calculating return value
 * @tparam RevGradFunT Type of tuple containing functors for reverse-mode
 *          gradients
 * @tparam FwdGradFunT Type of tuple containing functors for forward-mode
 *          gradients
 * @param input_input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param rev_grad_functors_tuple Tuple containing functors for reverse-mode
 *          gradients
 * @param fwd_grad_functors_tuple Tuple containing functors for forward-mode
 *          gradients
 * @return scalar or container of scalars with the values calculated by
 *          the value_functor functor and gradients calculated by the
 *          rev_grad_functors_tuple or fwd_grad_functors_tuple of functors,
 *          as appropriate
 */
template <typename ArgsTupleT, typename ValFunT, typename SharedArgsFunT,
          typename RevGradFunT, typename FwdGradFunT>
decltype(auto) function_gradients(ArgsTupleT&& input_args_tuple,
                                  ValFunT&& value_functor,
                                  SharedArgsFunT&& shared_args_functor,
                                  RevGradFunT&& rev_grad_functors_tuple,
                                  FwdGradFunT&& fwd_grad_functors_tuple) {
  return function_gradients_impl<true, return_type_t<ArgsTupleT>>(
      std::forward<ArgsTupleT>(input_args_tuple),
      std::forward<ValFunT>(value_functor),
      std::forward<SharedArgsFunT>(shared_args_functor),
      std::forward<RevGradFunT>(rev_grad_functors_tuple),
      std::forward<FwdGradFunT>(fwd_grad_functors_tuple));
}


/**
 * Specialisation for when function_gradients is called with only one tuple of
 * gradient functors, and a boolean flag indicating whether to use the functors
 * for both reverse and forward-mode.
 *
 * The forward-mode function_gradients_impl handles the different behaviours.
 *
 * @tparam ArgsTupleT Type of tuple containing all input arguments
 * @tparam ValFunT Type of functor for calculating return value
 * @tparam RevGradFunT Type of tuple containing functors for reverse-mode
 *          gradients
 * @tparam FwdGradFunT Type of tuple containing functors for forward-mode
 *          gradients
 * @param input_input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param grad_functors_tuple Tuple containing functors for gradients
 * @return scalar or container of scalars with the values calculated by
 *          the value_functor functor and gradients calculated by the
 *          grad_functors_tuple for reverse-mode and autodiff for
 *          forward-mode, as appropriate.
 */
template <bool FwdGradients, typename ArgsTupleT, typename ValFunT,
          typename SharedArgsFunT, typename GradFunT>
decltype(auto) function_gradients(ArgsTupleT&& input_args_tuple,
                                  ValFunT&& value_functor,
                                  SharedArgsFunT&& shared_args_functor,
                                  GradFunT&& grad_functors_tuple) {
    return function_gradients_impl<FwdGradients, return_type_t<ArgsTupleT>>(
        std::forward<ArgsTupleT>(input_args_tuple),
        std::forward<ValFunT>(value_functor),
        std::forward<SharedArgsFunT>(shared_args_functor),
        std::forward<GradFunT>(grad_functors_tuple),
        std::forward<GradFunT>(grad_functors_tuple));
}

/**
 * Specialisation for when function_gradients is called with only one tuple of
 * gradient functors, and no boolean flag indicating to use the functors for
 * both reverse and forward-mode.
 *
 * The gradient functors are assumed to be for reverse-mode, and that autodiff
 * should be used for forward-mode
 *
 * @tparam ArgsTupleT Type of tuple containing all input arguments
 * @tparam ValFunT Type of functor for calculating return value
 * @tparam RevGradFunT Type of tuple containing functors for reverse-mode
 *          gradients
 * @tparam FwdGradFunT Type of tuple containing functors for forward-mode
 *          gradients
 * @param input_input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param grad_functors_tuple Tuple containing functors for reverse-mode
 * @return scalar or container of scalars with the values calculated by
 *          the value_functor functor and gradients calculated by the
 *          grad_functors_tuple for reverse-mode and autodiff for
 *          forward-mode, as appropriate.
 */
template <typename ArgsTupleT, typename ValFunT, typename SharedArgsFunT,
          typename GradFunT>
decltype(auto) function_gradients(ArgsTupleT&& input_args_tuple,
                                  ValFunT&& value_functor,
                                  SharedArgsFunT&& shared_args_functor,
                                  GradFunT&& grad_functors_tuple) {
    return function_gradients_impl<false, return_type_t<ArgsTupleT>>(
        std::forward<ArgsTupleT>(input_args_tuple),
        std::forward<ValFunT>(value_functor),
        std::forward<SharedArgsFunT>(shared_args_functor),
        std::forward<GradFunT>(grad_functors_tuple),
        std::forward<GradFunT>(grad_functors_tuple));
}

}  // namespace math
}  // namespace stan

#endif

#ifndef STAN_MATH_FWD_FUNCTOR_FUNCTION_GRADIENTS_HPP
#define STAN_MATH_FWD_FUNCTOR_FUNCTION_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_tuple.hpp>
#include <stan/math/prim/functor/walk_tuples.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Helper function to create a variable for accumulating function gradients.
 * The function takes the return value as an input, so that the correct type
 * and dimensions can be deduced.
 * 
 * Specialisation for scalar inputs
 * 
 * @tparam T Type of return scalar
 * @param rtn_val Return value from function
 * @return Scalar of same type as return value, initialised to 0
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
constexpr T initialize_grad(T&& rtn_val) {
  return 0;
}

/**
 * Helper function to create a variable for accumulating function gradients.
 * 
 * Specialisation for Eigen inputs
 * 
 * @tparam T Type of return Eigen object
 * @param rtn_val Return value from function
 * @return Eigen object of same type as return value, initialised to all 0
 */
template <typename T, require_eigen_t<T>* = nullptr>
plain_type_t<T> initialize_grad(T&& rtn_eigen) {
  return plain_type_t<T>::Zero(rtn_eigen.rows(), rtn_eigen.cols());
}
}  // namespace internal

/**
 * Implementation functor for applying user-specified gradients for a function
 * with fvar<T> inputs.
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
 * @param input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param rev_grad_functors_tuple Tuple containing functors for reverse-mode
 *          gradients
 * @param fwd_grad_functors_tuple Tuple containing functors for forward-mode
 *          gradients
 * @return fvar<T> or container of fvar<T> with the values calculated by
 *          the value_functor functor and gradients by auto-differentiation
 */
template <bool FwdGradients, typename ScalarReturnT, typename ArgsTupleT,
          typename ValFunT, typename RevGradFunT, typename FwdGradFunT,
          require_st_fvar<ScalarReturnT>* = nullptr,
          require_not_t<std::integral_constant<bool, FwdGradients>>* = nullptr>
constexpr decltype(auto) function_gradients_impl(
    ArgsTupleT&& input_args_tuple, ValFunT&& value_functor,
    RevGradFunT&& rev_grad_functors_tuple, FwdGradFunT&& fwd_grad_functors_tuple) {
  return math::apply([&](auto&&... input_args) { return value_functor(input_args...); },
                    std::forward<ArgsTupleT>(input_args_tuple));
}

/**
 * Implementation functor for applying user-specified gradients for a function
 * with fvar<T> inputs.
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
 * @param input_args_tuple Tuple containing all input arguments
 * @param value_functor Functor for calculating return value
 * @param rev_grad_functors_tuple Tuple containing functors for reverse-mode
 *          gradients
 * @param fwd_grad_functors_tuple Tuple containing functors for forward-mode
 *          gradients
 * @return fvar<T> or container of fvar<T> with the values calculated by
 *          the value_functor functor and gradients by the fwd_grad_functors_tuple
 *          tuple of functors
 */
template <bool FwdGradients, typename ScalarReturnT, typename ArgsTupleT,
          typename ValFunT, typename RevGradFunT, typename FwdGradFunT,
          require_st_fvar<ScalarReturnT>* = nullptr,
          require_t<std::integral_constant<bool, FwdGradients>>* = nullptr>
decltype(auto) function_gradients_impl(
    ArgsTupleT&& input_args_tuple, ValFunT&& value_functor,
    RevGradFunT&& rev_grad_functors_tuple, FwdGradFunT&& fwd_grad_functors_tuple) {
  // Extract values from input arguments
  decltype(auto) values_tuple
      = map_tuple([&](auto&& input_arg) { return value_of(input_arg); },
                  std::forward<ArgsTupleT>(input_args_tuple));
  // Calculate return value
  decltype(auto) rtn_value
      = math::apply([&](auto&&... args) { return value_functor(args...); },
                    std::forward<decltype(values_tuple)>(values_tuple));

  // g++-4.9 has a bug with using decltype within a lambda for a captured
  // value. Create a tuple of unitialised values of the same type as the
  // function return that will be passed to the grad functor
  plain_type_t<decltype(rtn_value)> dummy_value;
  auto dummy_tuple = map_tuple([&](auto&& arg) { return dummy_value; },
                               std::forward<ArgsTupleT>(input_args_tuple));

  // Initialise variable for accumulating output gradients
  auto rtn_grad = internal::initialize_grad(std::forward<decltype(rtn_value)>(
    rtn_value));

  // Iterate over the tuple of input arguments and forward-mode gradient
  // functors, calculating the gradient at each and accumulating into the
  // rtn_grad variable
  walk_tuples(
      [&](auto&& grad_functor, auto&& input_arg, auto&& dummy) {
        using arg_t = decltype(input_arg);
        if (!is_constant_all<arg_t>::value) {
          rtn_grad += math::apply(
              [&](auto&&... args_values) {
                return grad_functor(rtn_value,
                         forward_as<promote_scalar_t<ScalarReturnT, arg_t>>(input_arg).d(),
                         args_values...);
              },
              values_tuple);
        }
      },
      std::forward<FwdGradFunT>(fwd_grad_functors_tuple),
      std::forward<ArgsTupleT>(input_args_tuple),
      std::forward<decltype(dummy_tuple)>(dummy_tuple));

  return to_fvar(rtn_value, rtn_grad);
}



}  // namespace math
}  // namespace stan
#endif

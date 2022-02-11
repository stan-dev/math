#ifndef STAN_MATH_PRIM_FUNCTOR_USER_GRADIENTS_HPP
#define STAN_MATH_PRIM_FUNCTOR_USER_GRADIENTS_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
/**
 * Specialisation for when user_gradients is called with all arithmetic
 * types.
 * 
 * @tparam ScalarT Return scalar type of function, used for SFINAE
 * @tparam ArgsTupleT Type of input tuple of arguments for function
 * @tparam ValFun Type of functor to calculate the return value
 * @tparam GradFunT Type of gradient functions tuple
 * @param args_tuple Tuple of input arguments to value functor
 * @param val_fun Functor for calculating return value
 * @param grad_fun_tuple Tuple of functors for calculating gradient
 * for each input
 * @return auto 
 */
template <typename ScalarT, typename ArgsTupleT,
          typename ValFun, typename GradFunT,
          require_arithmetic_t<ScalarT>* = nullptr>
auto user_gradients_impl(ArgsTupleT&& args_tuple, ValFun&& val_fun,
                         GradFunT&& grad_fun_tuple) {
  return math::apply([&](auto&&... args) {
    return val_fun(args...); }, std::forward<ArgsTupleT>(args_tuple)
  );
}

/**
 * Framework allowing users to provide gradient functions for their functions.
 * 
 * 
 * @tparam ArgsTupleT Type of input tuple of arguments for function
 * @tparam ValFun Type of functor to calculate the return value
 * @tparam GradFunT Type of gradient functions tuple
 * @param args_tuple Tuple of input arguments to value functor
 * @param val_fun Functor for calculating return value
 * @param grad_fun_tuple Tuple of functors for calculating gradient
 * for each input
 * @return auto 
 */
template <typename ArgsTupleT, typename ValFunT, typename GradFunTupleT>
auto user_gradients(const ArgsTupleT& args_tuple, const ValFunT& val_fun,
                    const GradFunTupleT& gradfun_tuple) {
  using scalar_rtn_t = scalar_type_t<return_type_t<ArgsTupleT>>;
  return user_gradients_impl<scalar_rtn_t>(args_tuple, val_fun, gradfun_tuple);
}

}  // namespace math
}  // namespace stan

#endif

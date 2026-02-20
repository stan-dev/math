#ifndef STAN_MATH_FWD_FUNCTOR_FINITE_DIFF_HPP
#define STAN_MATH_FWD_FUNCTOR_FINITE_DIFF_HPP

#include <stan/math/fwd/fun/tangent_of.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/zeroed_filtered_tuple.hpp>
#include <stan/math/prim/functor/finite_diff_gradient_auto.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/filter_types.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <tuple>

namespace stan {
namespace math {
namespace internal {
/**
 * Helper function for aggregating tangents if the respective input argument
 * was an fvar<T> type.
 *
 * Overload for when the input is not an fvar<T> and no tangents are needed.
 *
 * @tparam FuncTangent Type of tangent calculated by finite-differences
 * @tparam InputArg Type of the function input argument
 * @param[in] tangent Calculated tangent
 * @param[in] arg Input argument
 * @return `0.0`, because non-fvar inputs do not contribute to the tangent
 * aggregation.
 * @throw None.
 */
template <typename FuncTangent, typename InputArg,
          require_not_st_fvar<InputArg>* = nullptr>
inline constexpr double aggregate_tangent(const FuncTangent& tangent,
                                          const InputArg& arg) {
  return 0;
}

/**
 * Helper function for aggregating tangents if the respective input argument
 * was an fvar<T> type.
 *
 * Overload for when the input is an fvar<T> and its tangent needs to be
 * aggregated.
 *
 * @tparam FuncTangent Type of tangent calculated by finite-differences
 * @tparam InputArg Type of the function input argument
 * @param[in] tangent Calculated tangent
 * @param[in] arg Input argument
 * @return dot product of finite-difference tangent and input tangent.
 * @throw None.
 */
template <typename FuncTangent, typename InputArg,
          require_st_fvar<InputArg>* = nullptr>
inline auto aggregate_tangent(FuncTangent&& tangent, InputArg&& arg) {
  auto tangent_arr = as_array_or_scalar(std::forward<FuncTangent>(tangent));
  auto arg_tangent_arr
      = as_array_or_scalar(tangent_of(std::forward<InputArg>(arg)));
  if constexpr (is_stan_scalar_v<std::decay_t<decltype(tangent_arr)>>) {
    return tangent_arr * arg_tangent_arr;
  } else {
    using RetType = return_type_t<std::decay_t<decltype(tangent_arr(0))>,
                                  std::decay_t<decltype(arg_tangent_arr(0))>>;
    RetType rtn = 0;
    for (Eigen::Index i = 0; i < tangent_arr.size(); ++i) {
      rtn += tangent_arr(i) * arg_tangent_arr(i);
    }
    return rtn;
  }
}

}  // namespace internal

/**
 * Construct an fvar<T> where the tangent is calculated by finite-differencing.
 * Finite-differencing is only perfomed where the scalar type to be evaluated is
 * `fvar<T>.
 *
 * Higher-order inputs (i.e., fvar<var> & fvar<fvar<T>>) are also implicitly
 * supported through auto-diffing the finite-differencing process.
 *
 * @tparam F Type of functor for which fvar<T> support is needed
 * @tparam TArgs Template parameter pack of the types passed in the `operator()`
 *                of the functor type `F`. Must contain at least on type whose
 *                scalar type is `fvar<T>`
 *
 * Internal pattern used by this overload:
 * 1. Build `autodiff_args` with `filter_types<contains_autodiff>(args...)`.
 * 2. Build compact zero-initialized `grads` with
 *    `zeroed_filtered_tuple<contains_autodiff>(args...)`.
 * 3. Dispatch tuple-native finite differencing with:
 *    - `args_mask = args...` (original autodiff typing)
 *    - `args_work = value_of(args)...` (mutable finite-diff work tuple)
 * 4. Aggregate tangent contributions by zipping compact `grads` and
 *    `autodiff_args`.
 *
 * @param func Functor for which fvar<T> support is needed
 * @param args Parameter pack of arguments to be passed to functor.
 * @return Function value and tangent packed in the resulting `fvar`.
 * @throw Any exception thrown by `func`.
 */
template <typename F, typename... TArgs,
          require_any_st_fvar<TArgs...>* = nullptr>
inline auto finite_diff(const F& func, const TArgs&... args) {
  using FvarT = return_type_t<TArgs...>;
  using FvarInnerT = typename FvarT::Scalar;
  auto autodiff_args
      = stan::math::filter_types<contains_autodiff>(std::forward_as_tuple(args...));
  auto grads = stan::math::zeroed_filtered_tuple<contains_autodiff>(
      std::forward_as_tuple(args...));
  FvarInnerT rtn_value;
  finite_diff_gradient_auto(func, rtn_value, std::forward_as_tuple(args...),
                            std::forward_as_tuple(value_of(args)...), grads);
  FvarInnerT rtn_grad = 0;
  stan::math::for_each(
      [&rtn_grad](auto&& grad_i, auto&& arg_i) {
        rtn_grad += internal::aggregate_tangent(
            std::forward<decltype(grad_i)>(grad_i),
            std::forward<decltype(arg_i)>(arg_i));
      },
      grads, autodiff_args);

  return FvarT(rtn_value, rtn_grad);
}

/**
 * Construct an fvar<T> where the tangent is calculated by finite-differencing.
 * Finite-differencing is only perfomed where the scalar type to be evaluated is
 * `fvar<T>.
 *
 * This overload is used when no fvar<T> arguments are passed and simply
 * evaluates the functor with the provided arguments.
 *
 * @tparam F Type of functor
 * @tparam TArgs Template parameter pack of the types passed in the `operator()`
 *                of the functor type `F`. Must contain no type whose
 *                scalar type is `fvar<T>`
 * @param func Functor
 * @param args... Parameter pack of arguments to be passed to functor.
 * @return direct result of `func(args...)`.
 * @throw Any exception thrown by `func`.
 */
template <typename F, typename... TArgs,
          require_all_not_st_fvar<TArgs...>* = nullptr>
inline auto finite_diff(const F& func, const TArgs&... args) {
  return func(args...);
}

}  // namespace math
}  // namespace stan

#endif

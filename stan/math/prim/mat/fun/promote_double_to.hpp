#ifndef STAN_MATH_PRIM_MAT_FUN_PROMOTE_DOUBLE_TO_T_HPP
#define STAN_MATH_PRIM_MAT_FUN_PROMOTE_DOUBLE_TO_T_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/functor/apply.hpp>
#include <tuple>
#include <vector>
#include <iostream>

namespace stan {
namespace math {

template <typename T, typename... T_output>
auto promote_double_to_impl(std::tuple<T_output...> output) {
  return output;
}

template <typename T, int RowType, int ColType, typename... T_output,
          typename... T_inputs>
auto promote_double_to_impl(
    std::tuple<T_output...> output,
    const Eigen::Matrix<double, RowType, ColType>& leading_input,
    const T_inputs&... inputs);
template <typename T, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const std::vector<double>& leading_input,
                                   const T_inputs&... inputs);
template <typename T, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const double& leading_input,
                                   const T_inputs&... inputs);
template <typename T, typename R, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const R& leading_input,
                                   const T_inputs&... inputs);

template <typename T, typename R, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const R& leading_input,
                                   const T_inputs&... inputs) {
  return promote_double_to_impl<T>(output, inputs...);
}

template <typename T, int RowType, int ColType, typename... T_output,
          typename... T_inputs>
auto promote_double_to_impl(
    std::tuple<T_output...> output,
    const Eigen::Matrix<double, RowType, ColType>& leading_input,
    const T_inputs&... inputs) {
  Eigen::Matrix<T, RowType, ColType> promoted_leading_input(
      leading_input.rows(), leading_input.cols());
  for (int i = 0; i < leading_input.size(); ++i)
    promoted_leading_input(i) = leading_input(i);
  return promote_double_to_impl<T>(
      std::tuple_cat(output, std::make_tuple(promoted_leading_input)),
      inputs...);
}

template <typename T, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const std::vector<double>& leading_input,
                                   const T_inputs&... inputs) {
  std::vector<T> promoted_leading_input;
  promoted_leading_input.reserve(leading_input.size());
  for (int i = 0; i < leading_input.size(); ++i)
    promoted_leading_input.push_back(leading_input[i]);
  return promote_double_to_impl<T>(
      std::tuple_cat(output, std::make_tuple(promoted_leading_input)),
      inputs...);
}

template <typename T, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                                   const double& leading_input,
                                   const T_inputs&... inputs) {
  return promote_double_to_impl<T>(
      std::tuple_cat(output, std::make_tuple(T(leading_input))), inputs...);
}

/**
 * Cast the double elements of the given tuple to type T a return
 * them in another tuple. Drop the non-double input elements.
 *
 * @tparam T Type to cast double elements to
 * @tparam T_inputs Types in input tuple
 * @param input Input tuple
 * @return Tuple of double elements of input promoted to Ts
 */
template <typename T, typename... T_inputs>
auto promote_double_to(const std::tuple<T_inputs...>& input) {
  return apply([](auto ...lambda_args) {
      return promote_double_to_impl<T>(std::tuple<>(), lambda_args...);
    }, input);
}

}  // namespace math
}  // namespace stan
#endif

#ifndef STAN_MATH_PRIM_MAT_FUN_PROMOTE_DOUBLE_TO_T_HPP
#define STAN_MATH_PRIM_MAT_FUN_PROMOTE_DOUBLE_TO_T_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/functor/apply.hpp>
#include <tuple>
#include <vector>
#include <iostream>

namespace stan {
namespace math {

namespace internal {
template <typename T, typename... T_output>
auto promote_double_to_impl(std::tuple<T_output...> output) {
  return output;
}

template <typename T, int RowType, int ColType>
std::tuple<Eigen::Matrix<T, RowType, ColType> > promote_element_double_to(
    const Eigen::Matrix<double, RowType, ColType>& input) {
  Eigen::Matrix<T, RowType, ColType> promoted_input(input.rows(), input.cols());
  for (int i = 0; i < input.size(); ++i)
    promoted_input(i) = input(i);
  return std::make_tuple(promoted_input);
}

template <typename T>
std::tuple<std::vector<T> > promote_element_double_to(
    const std::vector<double>& input) {
  std::vector<T> promoted_input;
  promoted_input.reserve(input.size());
  for (size_t i = 0; i < input.size(); ++i)
    promoted_input.push_back(input[i]);
  return std::make_tuple(promoted_input);
}

template <typename T>
std::tuple<T> promote_element_double_to(const double& input) {
  return std::make_tuple(T(input));
}

template <typename T, typename R>
std::tuple<> promote_element_double_to(const R& input) {
  return std::tuple<>();
}

template <typename T, typename R, typename... T_output, typename... T_inputs>
auto promote_double_to_impl(std::tuple<T_output...> output,
                            const R& leading_input, const T_inputs&... inputs) {
  return promote_double_to_impl<T>(
      std::tuple_cat(output, promote_element_double_to<T>(leading_input)),
      inputs...);
}
}  // namespace internal

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
  return apply(
      [](auto... lambda_args) {
        return internal::promote_double_to_impl<T>(std::tuple<>(),
                                                   lambda_args...);
      },
      input);
}

}  // namespace math
}  // namespace stan
#endif

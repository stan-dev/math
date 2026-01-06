#ifndef STAN_MATH_REV_CORE_COLLECT_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_COLLECT_ADJOINTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta/is_all_arithmetic.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <type_traits>

namespace stan::math::internal {

inline constexpr bool ZeroOut = true;
/**
 * Collect the adjoints from the input and add them to the output.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_var_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                          input_i.size());
          output_map.array() += input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i.adj();
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) collect_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

/**
 * Collects the adjoints from the input and adds them to the output.
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_arithmetic_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output&& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<double, -1, 1>> input_map(input_i.data(),
                                                             input_i.size());
          output_map.array() += input_map.array();
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.array();
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i;
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) collect_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

/**
 * Used in reverse pass to collect adjoints to the output
 * @tparam Output A tuple or type where all scalar types are `var` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param ret The vari object containing the adjoint to be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void collect_adjoints(Output&& output, const vari* ret, Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    static_assert(sizeof(std::decay_t<Output>*) == 0,
                  "INTERNAL ERROR:(laplace_marginal_lpdf) "
                  "Accumulate Adjoints called on a tuple, but tuples cannot be "
                  "on the reverse mode stack! "
                  "This is an internal error, please report it: "
                  "https://github.com/stan-dev/math/issues");
  } else if constexpr (is_std_vector_v<Output>) {
    if constexpr (!is_var_v<value_type_t<Output>>) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], ret, input[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<var, -1, 1>> output_map(output.data(),
                                                       output.size());
      Eigen::Map<const Eigen::Matrix<double, -1, 1>> input_map(input.data(),
                                                               input.size());
      output_map.array().adj() += ret->adj_ * input_map.array();
    }
  } else if constexpr (is_eigen_v<Output>) {
    output.adj().array() += ret->adj_ * input.array();
  } else if constexpr (is_var_v<Output>) {
    output.adj() += ret->adj_ * input;
  }
}

}  // namespace stan::math::internal

#endif

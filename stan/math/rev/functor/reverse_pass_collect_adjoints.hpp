#ifndef STAN_MATH_REV_FUNCTOR_REVERSE_PASS_COLLECT_ADJOINTS_HPP
#define STAN_MATH_REV_FUNCTOR_REVERSE_PASS_COLLECT_ADJOINTS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/collect_adjoints.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <cstddef>
#include <utility>

namespace stan::math::internal {

/**
 * Collects adjoints from a tuple or std::vector of tuples
 * @tparam Output A tuple or std::vector of tuples where all scalar types are
 * `var` types
 * @tparam Input A tuple or std::vector of tuples where all scalar types are
 * `arithmetic` types
 * @param ret The vari object containing the adjoint to be added
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void reverse_pass_collect_adjoints(var ret, Output&& output,
                                          Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    stan::math::for_each(
        [ret](auto&& inner_arg, auto&& inner_input) mutable {
          reverse_pass_collect_adjoints(
              ret, std::forward<decltype(inner_arg)>(inner_arg),
              std::forward<decltype(inner_input)>(inner_input));
        },
        std::forward<Output>(output), std::forward<Input>(input));
  } else if constexpr (is_std_vector_containing_tuple_v<Output>) {
    for (std::size_t i = 0; i < output.size(); ++i) {
      reverse_pass_collect_adjoints(ret, output[i], input[i]);
    }
  } else {
    reverse_pass_callback(
        [vi = ret.vi_, arg_arena = to_arena(std::forward<Output>(output)),
         input_arena = to_arena(std::forward<Input>(input))]() mutable {
          collect_adjoints(arg_arena, vi, input_arena);
        });
  }
}

}  // namespace stan::math::internal

#endif

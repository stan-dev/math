#ifndef STAN_MATH_REV_FUN_BUILD_MATRIX
#define STAN_MATH_REV_FUN_BUILD_MATRIX

#include <stan/math/prim/fun/build_matrix.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <stan/math/prim/functor/map_if.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/from_var_value.hpp>
#include <stan/math/rev/fun/to_arena.hpp>

namespace stan::math {

namespace internal {

/**
 * Deep copy the var-containing arguments in a tuple onto the nested autodiff
 * stack while forwarding arithmetic-only arguments by reference.
 *
 * @tparam ArgsTuple tuple of arguments
 * @param[in] args_tup tuple of arguments to copy
 * @return tuple with var-containing arguments deep copied
 */
template <typename ArgsTuple, require_tuple_t<ArgsTuple>* = nullptr>
inline auto deep_copy_vars_tuple(ArgsTuple&& args_tup) {
  return map_if<is_any_var_scalar>(
      [](auto&& arg) {
        return deep_copy_vars(std::forward<decltype(arg)>(arg));
      },
      std::forward<ArgsTuple>(args_tup));
}

/**
 * Build a tuple of value-only references suitable for primal evaluation.
 *
 * Each entry is converted with `value_of` and wrapped in `to_ref` so Eigen
 * expressions are evaluated only once and reused across the matrix fill.
 *
 * @tparam ArgsTuple tuple of arguments
 * @param[in] args_tup tuple of arguments to convert
 * @return tuple of value-only references
 */
template <typename ArgsTuple, require_tuple_t<ArgsTuple>* = nullptr>
inline auto make_value_tuple(ArgsTuple&& args_tup) {
  return apply(
      [](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      std::forward<ArgsTuple>(args_tup));
}

/**
 * Construct a matrix by evaluating `f(i, j, args...)` at each entry.
 *
 * @tparam F matrix element functor
 * @tparam ArgsTuple tuple of arguments passed to `f`
 * @param[in] rows number of rows to evaluate
 * @param[in] cols number of columns to evaluate
 * @param[in] f callable used to generate each matrix element
 * @param[in] args_tup tuple of arguments forwarded to `f`
 * @return matrix of evaluated entries
 */
template <typename F, typename ArgsTuple, require_tuple_t<ArgsTuple>* = nullptr>
inline auto fill_matrix(int rows, int cols, F&& f, ArgsTuple&& args_tup) {
  using scalar_t = return_type_t<ArgsTuple>;
  return apply(
      [rows, cols, &f](auto&&... args) {
        Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> out(rows, cols);
        for (int i = 0; i < rows; ++i) {
          for (int j = 0; j < cols; ++j) {
            out(i, j) = f(i, j, args...);
          }
        }
        return out;
      },
      std::forward<ArgsTuple>(args_tup));
}

}  // namespace internal

/**
 * Construct a reverse-mode autodiff matrix by evaluating `f(i, j, args...)`
 * for each element.
 *
 * The forward pass computes only the primal matrix values. The reverse pass
 * rebuilds the matrix once on a nested autodiff tape, seeds the nested matrix
 * adjoints with the outer matrix adjoints, and performs a single nested
 * reverse sweep to accumulate adjoints back to the input arguments.
 *
 * @tparam F matrix element functor
 * @tparam Args argument pack to `f`
 * @param[in] rows number of matrix rows
 * @param[in] cols number of matrix columns
 * @param[in] f callable used to generate each matrix element
 * @param[in] args arguments forwarded to `f`
 * @return matrix of vars whose reverse pass applies the matrix pullback once
 */
template <typename F, typename... Args, require_any_st_var<Args...>* = nullptr>
inline auto build_matrix(int rows, int cols, F&& f, Args&&... args) {
  auto arena_args_tuple = to_arena(std::make_tuple(eval(args)...));
  auto args_vals_tuple = internal::make_value_tuple(arena_args_tuple);
  const size_t num_vars = count_vars(args...);

  arena_t<Eigen::MatrixXd> vals
      = internal::fill_matrix(rows, cols, f, args_vals_tuple);

  auto out = make_callback_var(
      std::move(vals),
      [rows, cols, num_vars, arena_args_tuple,
       fun = std::decay_t<F>(std::forward<F>(f))](auto& vi) mutable {
        if (vi.size() == 0 || num_vars == 0) {
          return;
        }

        nested_rev_autodiff nested;
        auto args_copy = internal::deep_copy_vars_tuple(arena_args_tuple);
        auto args_vars_copy = internal::filter_var_scalar_types(args_copy);
        auto args_vars_refs
            = internal::filter_var_scalar_types(arena_args_tuple);
        auto nested_m = internal::fill_matrix(rows, cols, fun, args_copy);
        nested_m.adj() = vi.adj();
        grad();
        iter_tuple_nested(
            [](auto&& output_i, auto&& input_i) {
              using output_t = std::decay_t<decltype(output_i)>;
              if constexpr (is_std_vector_v<output_t>) {
                Eigen::Map<Eigen::Matrix<var, -1, 1>> output_map(
                    output_i.data(), output_i.size());
                Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                                input_i.size());
                output_map.adj().array() += input_map.adj().array();
              } else if constexpr (is_eigen_v<output_t>) {
                output_i.adj().array() += input_i.adj().array();
              } else if constexpr (is_var<output_t>::value) {
                output_i.adj() += input_i.adj();
              } else {
                static_assert(
                    sizeof(std::decay_t<output_t>*) == 0,
                    "INTERNAL ERROR: build_matrix encountered an unsupported "
                    "type. This is an internal error, please report it: "
                    "https://github.com/stan-dev/math/issues");
              }
            },
            args_vars_refs, args_vars_copy);
      });
  return from_var_value(out);
}
}  // namespace stan::math

#endif

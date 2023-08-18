#ifndef TEST_UNIT_MATH_EXPR_TESTS_HPP
#define TEST_UNIT_MATH_EXPR_TESTS_HPP

#include <stan/math/mix.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/expressions/expression_test_helpers.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace stan {
namespace test {
namespace internal {

template <typename T>
struct char_scalar_type;

template <>
struct char_scalar_type<double> {
  static constexpr const char* scalar{"double"};
};
constexpr const char* char_scalar_type<double>::scalar;

template <>
struct char_scalar_type<stan::math::var> {
  static constexpr const char* scalar{"var"};
};
constexpr const char* char_scalar_type<stan::math::var>::scalar;

template <>
struct char_scalar_type<stan::math::fvar<double>> {
  static constexpr const char* scalar{"fvar<double>"};
};
constexpr const char* char_scalar_type<stan::math::fvar<double>>::scalar;

template <>
struct char_scalar_type<std::complex<double>> {
  static constexpr const char* scalar{"std::complex<double>"};
};
constexpr const char* char_scalar_type<std::complex<double>>::scalar;

template <>
struct char_scalar_type<std::complex<stan::math::var>> {
  static constexpr const char* scalar{"std::complex<var>"};
};
constexpr const char* char_scalar_type<std::complex<stan::math::var>>::scalar;

template <>
struct char_scalar_type<std::complex<stan::math::fvar<double>>> {
  static constexpr const char* scalar{"std::complex<fvar<double>>"};
};
constexpr const char*
    char_scalar_type<std::complex<stan::math::fvar<double>>>::scalar;

/**
 * Check that the evaluations were less than the size of the input size.
 *
 * For the expression tests, this function is used to check that each input
 * is not accessed more times than the number of elements available in the
 * Eigen matrix inputs. Each element of the two input arrays correspond to
 * an argument passed to the function that is being checked with the expression
 * tests. All elements of an array that do not correspond to an Eigen matrix
 * input will be set to zero and ignored.
 *
 * @param arg_evals How many times an eigen matrix argument was evaluated.
 *  For non Eigen matrix types the value will be 0.
 * @param size_of_arg The size of a Eigen matrix type. For non eigen types
 *  the associated array will be zero.
 */
template <typename ScalarType, std::size_t N>
void expect_all_used_only_once(std::array<int, N>& arg_evals,
                               std::array<int, N>& size_of_arg) {
  for (int i = 0; i < N; ++i) {
    EXPECT_LE(arg_evals[i], size_of_arg[i])
        << "(" << char_scalar_type<ScalarType>::scalar << ")"
        << " argument " << std::to_string(i) << " was evaluated "
        << std::to_string(arg_evals[i])
        << " times but should"
           " be evaluated no more than its size ("
        << std::to_string(size_of_arg[i])
        << "). Before accessing an input matrix in a function use "
           "\n`auto&& {input_name}_ref = stan::math::to_ref({input_name})`\n"
           " to evaluate the input matrix once."
           " Then use `{input_name}_ref` inside of the function instead"
           " of the original input.";
  }
}

/**
 * For eigen types, construct a `CounterOp` expression to keep track of access.
 *
 * `stan::test::counterOp<stan::scalar_type_t<EigMat>>` is a Eigen UnaryOp
 * class which keeps track of a pointer to an integer. Everytime the CounterOp's
 * `operator()` is called the integer gets ticked up by one value. We use
 * `CounterOp` as a way to test if an expression is being called more times than
 *  the number of elements in the matrix. If the number of calls exceeds the
 *  * size of the matrix then we have run the expression too many times.
 *
 * @tparam EigMat A type derived from `Eigen::EigenBase`
 * @param count An integer used to keep track of the number of expression
 *  evaluations.
 * @param arg A type derived from `Eigen::EigenBase`
 */
template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto make_expr(int& count, EigMat&& arg) {
  return arg.unaryExpr(
      stan::test::counterOp<stan::scalar_type_t<EigMat>>(&count));
}

/**
 * For non-eigen types this function just returns `arg`
 * @tparam Any type not derived from `Eigen::EigenBase`
 */
template <typename T, stan::require_not_eigen_t<T>* = nullptr>
auto make_expr(int& /* count */, T&& arg) {
  return arg;
}

/**
 * Generate a tuple from input args where Eigen types are wrapped by `counterOp`
 * @tparam N The number of arguments
 * @tparam Args Parameter pack of arguments passed to function to test.
 * @param expr_evals Integer array used to keep track of accesses to each
 * argument
 * @param args Arguments to pass to function. Eigen arguments will be wrapped in
 *  `counterOp`.
 */
template <std::size_t N, typename... Args>
auto make_expr_args(std::array<int, N>& expr_evals, Args&&... args) {
  return stan::math::index_apply<N>([&expr_evals, &args...](auto... Is) {
    return std::make_tuple(make_expr(expr_evals[Is], args)...);
  });
}

/**
 * For non-eigen types return 0.
 */
template <typename T, stan::require_not_eigen_t<T>* = nullptr>
inline constexpr int eigen_size(T&& x) {
  return 0;
}

/**
 * Get the size of an eigen type.
 * @tparam EigMat A type derived from `Eigen::EigenBase`
 */
template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
inline int eigen_size(EigMat&& x) {
  return x.size();
}

/**
 * Return an integer array with counts of the eigen matrices in a parameter pack
 * @tparam Args A parameter pack of argument types
 * @param args A parameter pack of arguments.
 * @return An integer array where each element corresponds the size of
 *  type in the associated `args` position. If the associated argument is not
 *  an Eigen type then the value is be zero.
 */
template <typename... Args>
std::array<int, sizeof...(Args)> eigen_arg_sizes(Args&&... args) {
  return std::array<int, sizeof...(Args)>{eigen_size(args)...};
}

/**
 * No-op for when none of the arguments is an Eigen type
 */
template <typename ScalarType, typename F, typename... Args,
          require_all_not_eigen_t<Args...>* = nullptr>
void check_expr_test(F&& f, Args&&... args) {}

/**
 * Check whether any Eigen inputs are executed too many times.
 *
 * This function checks whether a function is safe for passing in
 *  Eigen expressions by evaluating whether a passed in expression would be
 *  evaluated more than the number of elements available in the Eigen
 *  expression's underlying matrix. Each Eigen object passed to this function
 *  is wrapped in a `counterOp` unary expression which keeps track of the
 *  number of times it's `operator()` is called. If the `counterOp` is called
 *  more times than the number of elements in the underlying matrix, then
 *  excess computation has been done and so the test will fail.
 *
 * Note that this does not catch all possible failures. For instance a
 * function that only adds up the top half of a vector twice would pass this
 * check, since the number of accesses would be smaller or equal to the total
 * size of the vector. But that would fail the spirit of the test, which is to
 * make sure we are not needlessly evaluating expressions multiple times.
 *
 *
 * @tparam ScalarType Passed by user explicitly and decides the scalar type
 *  to promote the arguments in `Args` to.
 * @tparam F a functor with an `operator(Args&&...)`
 * @tparam Args Parameter pack with at least one Eigen type.
 * @param f functor whose `operator()` will be called.
 * @param args pack of arguments to pass to the functor.
 */
template <typename ScalarType, typename F, typename... Args,
          require_any_eigen_t<Args...>* = nullptr>
void check_expr_test(F&& f, Args&&... args) {
  std::array<int, sizeof...(args)> expr_eval_counts;
  for (int i = 0; i < sizeof...(args); ++i) {
    expr_eval_counts[i] = 0;
  }
  std::array<int, sizeof...(args)> size_of_eigen_args
      = eigen_arg_sizes(args...);
  // returns tuple of unary count expressions
  auto promoted_args = std::make_tuple(
      stan::math::eval(stan::math::promote_scalar<ScalarType>(args))...);
  auto expr_args = stan::math::apply(
      [&expr_eval_counts](auto&&... args) mutable {
        return make_expr_args(expr_eval_counts, args...);
      },
      promoted_args);
  auto return_val = stan::math::eval(stan::math::apply(
      [&f](auto&&... args) { return f(std::forward<decltype(args)>(args)...); },
      expr_args));
  expect_all_used_only_once<ScalarType>(expr_eval_counts, size_of_eigen_args);
  if (stan::is_var<ScalarType>::value) {
    stan::math::recover_memory();
  }
}

}  // namespace internal
/**
 * Check whether any Eigen inputs are executed too many times.
 *
 * @tparam F a functor with an `operator(Args&&...)`
 * @tparam Args Parameter pack with at least one Eigen type.
 * @param f functor whose `operator()` will be called.
 * @param args pack of arguments to pass to the functor.
 */
template <typename F, typename... Args,
          require_all_st_stan_scalar<Args...>* = nullptr,
          require_all_not_st_complex<Args...>* = nullptr>
void check_expr_test(F&& f, Args&&... args) {
  try {
    stan::test::internal::check_expr_test<double>(f, args...);
    try {
      stan::test::internal::check_expr_test<stan::math::var>(f, args...);
      stan::math::recover_memory();
    } catch (const std::exception& e) {
      stan::math::recover_memory();
    }
    stan::test::internal::check_expr_test<stan::math::fvar<double>>(f, args...);
  } catch (const std::exception& e) {
  }
}

template <typename F, typename... Args,
          require_any_st_complex<Args...>* = nullptr>
void check_expr_test(F&& f, Args&&... args) {
  try {
    stan::test::internal::check_expr_test<std::complex<double>>(f, args...);
    try {
      stan::test::internal::check_expr_test<std::complex<stan::math::var>>(
          f, args...);
      stan::math::recover_memory();
    } catch (const std::exception& e) {
      stan::math::recover_memory();
    }
    stan::test::internal::check_expr_test<
        std::complex<stan::math::fvar<double>>>(f, args...);
  } catch (const std::exception& e) {
  }
}

}  // namespace test
}  // namespace stan

#endif

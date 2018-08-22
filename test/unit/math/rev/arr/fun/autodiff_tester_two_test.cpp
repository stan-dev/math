#include <stan/math/rev/arr.hpp>
#include <stan/math/prim/scal/fun/call_all_argument_combos.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

template <typename F, typename... Targs>
void test_gradients(F f, double dx, const Targs&... args) {
  stan::math::start_nested();

  auto fd = [&f](auto... args) { return f(stan::math::value_of(args)...); };
  auto input = stan::math::variable_adapter_factory<stan::math::var>(args...);
  auto output = stan::math::variable_adapter_factory<stan::math::var>(
      std::apply(f, input.args_));

  if (input.size() == 0 || output.size() == 0) {
    return;
  }

  Eigen::MatrixXd fd_jac(output.size(), input.size());
  for (size_t j = 0; j < input.size(); ++j) {
    stan::math::variable_adapter fd_input = input;
    fd_input(j) += dx;
    auto output1 = stan::math::variable_adapter_factory<double>(
        std::apply(fd, fd_input.args_));
    fd_input(j) -= 2 * dx;
    auto output2 = stan::math::variable_adapter_factory<double>(
        std::apply(fd, fd_input.args_));
    for (size_t i = 0; i < output.size(); ++i) {
      fd_jac(i, j) = (output1(i) - output2(i)) / (2.0 * dx);
    }
  }

  Eigen::MatrixXd var_jac(output.size(), input.size());
  for (size_t i = 0; i < output.size(); ++i) {
    output(i).grad();

    for (size_t j = 0; j < input.size(); ++j) {
      var_jac(i, j) = input(j).adj();
    }

    stan::math::set_zero_all_adjoints();
  }

  for (int i = 0; i < fd_jac.rows(); ++i) {
    for (int j = 0; j < fd_jac.cols(); ++j) {
      EXPECT_NEAR(fd_jac(i, j), var_jac(i, j), dx * dx * 1e2);
    }
  }

  stan::math::recover_memory_nested();
}

template <typename F, typename... Targs>
void autodiff_tester_two(F f, const Targs&... args) {
  stan::math::call_all_argument_combos(
      [&f](auto... args) {
        test_gradients(f, 1e-4, args...);
        return 0;
      },
      std::tuple_cat(
          args, stan::math::promote_double_to_T<stan::math::var>(args))...);
}

TEST(MathFunctions, autodiff_tester_two_basic_function) {
  auto func = [](auto a, auto b) { return a + b; };

  autodiff_tester_two(func, std::make_tuple(1, 1.0), std::make_tuple(2, 3.0));
}

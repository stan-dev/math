#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

namespace apply_test {
struct func {
  template <typename T>
  T operator()(T t) {
    return t;
  }
};
}  // namespace apply_test

TEST(MathFunctions, apply_basic_empty) {
  std::tuple<> x;

  auto y = stan::math::apply([]() { return static_cast<double>(1.7); }, x);

  EXPECT_EQ(1.7, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_basic_double) {
  std::tuple<double> x = std::make_tuple(1.0);

  auto y = stan::math::apply(apply_test::func{}, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_basic_int) {
  std::tuple<int> x = std::make_tuple(1);

  auto y = stan::math::apply(apply_test::func{}, x);

  EXPECT_EQ(1, y);
  EXPECT_TRUE((std::is_same<int, decltype(y)>::value));
}

TEST(MathFunctions, apply_const_double) {
  const std::tuple<const double> x = std::make_tuple(1.0);

  auto y = stan::math::apply(apply_test::func{}, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_double) {
  auto y = stan::math::apply(apply_test::func{}, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function) {
  auto y = stan::math::apply([](auto x) { return x; }, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_reference) {
  auto y = stan::math::apply([](auto& x) { return x; }, std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_const_reference) {
  auto y = stan::math::apply([](const auto& x) { return x; },
                             std::make_tuple(1.0));

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_return_reference) {
  std::tuple<double> x = std::make_tuple(1.0);

  decltype(auto) y = stan::math::apply(
      [](auto& x) -> auto& { return x; }, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<double&, decltype(y)>::value));
}

TEST(MathFunctions, apply_temporary_function_return_const_reference) {
  std::tuple<double> x = std::make_tuple(1.0);

  decltype(auto) y = stan::math::apply(
      [](const auto& x) -> const auto& { return x; }, x);

  EXPECT_EQ(1.0, y);
  EXPECT_TRUE((std::is_same<const double&, decltype(y)>::value));
}

TEST(MathFunctions, apply_two_args) {
  std::tuple<double, int> x = std::make_tuple(1.5, 1);

  auto y = stan::math::apply([](double x0, int x1) { return x0 + x1; }, x);

  EXPECT_EQ(2.5, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

TEST(MathFunctions, apply_three_args) {
  std::tuple<double, double, int> x = std::make_tuple(1.5, 0.75, 1);

  auto y = stan::math::apply(
      [](double x0, double x1, int x2) { return x0 + x1 + x2; }, x);

  EXPECT_EQ(3.25, y);
  EXPECT_TRUE((std::is_same<double, decltype(y)>::value));
}

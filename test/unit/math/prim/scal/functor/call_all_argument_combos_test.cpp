#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>

TEST(MathFunctions, call_all_argument_combos_0_arg) {
  auto output = stan::math::call_all_argument_combos([]() { return 0; });

  EXPECT_EQ(1, std::tuple_size<decltype(output)>::value);
  EXPECT_EQ(0, std::get<0>(output));
}

TEST(MathFunctions, call_all_argument_combos_2_arg) {
  auto output = stan::math::call_all_argument_combos(
      [](auto a, auto b) { return a + b; }, std::make_tuple(1, 2),
      std::make_tuple(3, 4));

  EXPECT_EQ(4, std::tuple_size<decltype(output)>::value);
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<0>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<1>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<2>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<3>(output))>::value));

  EXPECT_EQ(4, std::get<0>(output));
  EXPECT_EQ(5, std::get<1>(output));
  EXPECT_EQ(5, std::get<2>(output));
  EXPECT_EQ(6, std::get<3>(output));
}

TEST(MathFunctions, call_all_argument_combos_3_arg) {
  auto output = stan::math::call_all_argument_combos(
      [](auto a, auto b, auto c) { return a + b + c; }, std::make_tuple(1, 2),
      std::make_tuple(3, 4), std::make_tuple(5, 6));

  EXPECT_EQ(8, std::tuple_size<decltype(output)>::value);
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<0>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<1>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<2>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<3>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<4>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<5>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<6>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<7>(output))>::value));

  EXPECT_EQ(9, std::get<0>(output));
  EXPECT_EQ(10, std::get<1>(output));
  EXPECT_EQ(10, std::get<2>(output));
  EXPECT_EQ(11, std::get<3>(output));
  EXPECT_EQ(10, std::get<4>(output));
  EXPECT_EQ(11, std::get<5>(output));
  EXPECT_EQ(11, std::get<6>(output));
  EXPECT_EQ(12, std::get<7>(output));
}

TEST(MathFunctions, call_all_argument_combos_2_arg_types) {
  auto output = stan::math::call_all_argument_combos(
      [](auto a, auto b) { return a + b; }, std::make_tuple(1.1, 2),
      std::make_tuple(3, 4.7));

  EXPECT_EQ(4, std::tuple_size<decltype(output)>::value);
  EXPECT_TRUE((std::is_same<double&, decltype(std::get<0>(output))>::value));
  EXPECT_TRUE((std::is_same<double&, decltype(std::get<1>(output))>::value));
  EXPECT_TRUE((std::is_same<int&, decltype(std::get<2>(output))>::value));
  EXPECT_TRUE((std::is_same<double&, decltype(std::get<3>(output))>::value));

  EXPECT_FLOAT_EQ(4.1, std::get<0>(output));
  EXPECT_FLOAT_EQ(5.8, std::get<1>(output));
  EXPECT_EQ(5, std::get<2>(output));
  EXPECT_FLOAT_EQ(6.7, std::get<3>(output));
}

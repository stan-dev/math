#include <stan/math/prim.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

TEST(MathFunctions, iter_tuple_nested_empty) {
  auto x = 1;
  stan::math::iter_tuple_nested([&x](auto&& args) { return x++; },
                                std::make_tuple());
  // This should never call the lambda, so x should not change.
  EXPECT_EQ(x, 1);
}

TEST(MathFunctions, iter_tuple_nested_basic) {
  auto output = std::make_tuple(1, 2, 3);
  stan::math::iter_tuple_nested([](auto&& arg1, auto&& arg2) { arg1 += arg2; },
                                output, std::make_tuple(1, 2, 3));
  EXPECT_EQ(std::get<0>(output), 2);
  EXPECT_EQ(std::get<1>(output), 4);
  EXPECT_EQ(std::get<2>(output), 6);
}

TEST(MathFunctions, iter_tuple_nested_deep_tuple) {
  using inner_vec_t = std::vector<std::tuple<int, int>>;
  inner_vec_t inner_val{{1, 2}, {2, 3}};
  using inner_t = std::tuple<inner_vec_t, inner_vec_t>;
  auto output = std::make_tuple(1, inner_t{inner_val, inner_val}, 3);
  auto input = output;
  stan::math::iter_tuple_nested(
      [](auto&& arg1, auto&& arg2) { return arg1 += arg2; }, output, input);
  EXPECT_EQ(std::get<0>(output), 2);

  auto&& inner_output = std::get<0>(std::get<1>(output));
  auto&& inner_val_i = inner_output[0];
  EXPECT_EQ(std::get<0>(inner_val_i), 2);
  EXPECT_EQ(std::get<1>(inner_val_i), 4);
  inner_val_i = inner_output[1];
  EXPECT_EQ(std::get<0>(inner_val_i), 4);
  EXPECT_EQ(std::get<1>(inner_val_i), 6);

  inner_output = std::get<1>(std::get<1>(output));
  inner_val_i = inner_output[0];
  EXPECT_EQ(std::get<0>(inner_val_i), 2);
  EXPECT_EQ(std::get<1>(inner_val_i), 4);
  inner_val_i = inner_output[1];
  EXPECT_EQ(std::get<0>(inner_val_i), 4);
  EXPECT_EQ(std::get<1>(inner_val_i), 6);

  EXPECT_EQ(std::get<2>(output), 6);
}

}  // namespace

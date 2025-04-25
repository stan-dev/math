#include <stan/math/prim.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

TEST(MathFunctions, filter_map_empty) {
  auto args = stan::math::filter_map<stan::test::always_true>(
      [](auto&& arg) { return arg; }, std::make_tuple());
  EXPECT_EQ(std::tuple_size_v<decltype(args)>, 0);
}

TEST(MathFunctions, filter_map_all_true_filter) {
  auto args = stan::math::filter_map<stan::test::always_true>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::make_tuple(1, 2, 3));
  static_assert(std::tuple_size_v<decltype(args)> == 3,
                "tuple size should be 3!");
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(args)>, int>));
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, decltype(args)>, int>));
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<2, decltype(args)>, int>));
}
TEST(MathFunctions, filter_map_all_false) {
  auto args = stan::math::filter_map<stan::test::always_false>(
      [](auto&& arg) { return arg; }, std::make_tuple(1, 2, 3));
  static_assert(std::tuple_size_v<decltype(args)> == 0,
                "tuple size should be 0!");
}

TEST(MathFunctions, filter_map_first_true) {
  auto args = stan::math::filter_map<std::is_floating_point>(
      [](auto&& arg) { return arg; },
      std::make_tuple(1.0, 2.0, 3, 4.0, 5, 6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 4,
                "tuple size should be 4!");
}

TEST(MathFunctions, filter_map_first_false) {
  auto args = stan::math::filter_map<std::is_floating_point>(
      [](auto&& arg) { return arg; }, std::make_tuple(1, 2, 3, 4.0, 5, 6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 2,
                "tuple size should be 2!");
}

TEST(MathFunctions, filter_map_inner_tuple_rvalue) {
  auto args = stan::math::filter_map<stan::test::contains_floating_point>(
      [](auto&& arg) { return arg; },
      std::forward_as_tuple(1, 2, 3, 4.0, 5, std::forward_as_tuple(1.0, 2, 3.0),
                            6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 3,
                "tuple size should be 3!");
  using inner_tuple_t = std::tuple_element_t<1, decltype(args)>;
  static_assert(!stan::test::is_ref_element_v<0, decltype(args)>,
                "0th should not be a reference!");
  static_assert(!stan::test::is_ref_element_v<1, decltype(args)>,
                "1st should not be a reference!");
  static_assert(!stan::test::is_ref_element_v<0, inner_tuple_t>,
                "tuple<1><0> should not be a reference!");
  static_assert(!stan::test::is_ref_element_v<1, inner_tuple_t>,
                "tuple<1><1> should not be a reference!");
  static_assert(!stan::test::is_ref_element_v<2, decltype(args)>,
                "2nd should not be a reference!");
}

TEST(MathFunctions, filter_map_inner_tuple_mix_value_type) {
  using stan::test::is_lvalue_ref_element_v;
  using stan::test::is_same_tuple_element_v;
  double a = 1.0;
  std::vector<int> b{1, 2, 3};
  auto args = stan::math::filter_map<stan::test::is_fp_or_std_vector>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::forward_as_tuple(
          1, 2.0, b, 4, 5, a, std::vector<int>{0, 8, 7},
          std::forward_as_tuple(1, std::vector<int>{3, 2, 1}, b), b, 6.0));
  using args_t = decltype(args);
  static_assert(is_same_tuple_element_v<0, args_t, double>,
                "0th should be a double!");
  static_assert(!stan::test::is_ref_element_v<0, args_t>,
                "0th should not be a reference!");
  static_assert(is_lvalue_ref_element_v<1, args_t>, "1st should be an lvalue!");
  static_assert(is_same_tuple_element_v<1, args_t, std::vector<int>>,
                "1st should be an std::vector<int>!");
  static_assert(is_lvalue_ref_element_v<2, args_t>, "2nd should be an lvalue!");
  static_assert(is_same_tuple_element_v<2, args_t, double>,
                "2nd should be a double!");
  static_assert(!stan::test::is_ref_element_v<3, args_t>,
                "3rd should not be a reference!");
  static_assert(is_same_tuple_element_v<3, args_t, std::vector<int>>,
                "3rd should be an std::vector<int>!");
  using inner_tuple_t = std::tuple_element_t<4, args_t>;
  static_assert(!stan::test::is_ref_element_v<0, inner_tuple_t>,
                "tuple<4><0> should not be a reference!");
  static_assert(is_same_tuple_element_v<0, inner_tuple_t, std::vector<int>>,
                "tuple<4><0> should be an std::vector<int>!");
  static_assert(is_lvalue_ref_element_v<1, inner_tuple_t>,
                "tuple<4><1> should be an lvalue!");
  static_assert(is_same_tuple_element_v<1, inner_tuple_t, std::vector<int>>,
                "tuple<4><1> should be an std::vector<int>!");
  static_assert(is_lvalue_ref_element_v<5, args_t>, "5th should be an lvalue!");
  static_assert(is_same_tuple_element_v<5, args_t, std::vector<int>>,
                "5th should be an std::vector<int>!");
  static_assert(!stan::test::is_ref_element_v<6, args_t>,
                "6th should not be a reference!");
  static_assert(is_same_tuple_element_v<6, args_t, double>,
                "6th should be a double!");
}

}  // namespace

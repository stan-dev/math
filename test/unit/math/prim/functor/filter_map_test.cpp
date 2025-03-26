#include <stan/math/prim.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {
template <typename T>
struct always_true {
  static constexpr bool value = true;
};
template <typename T>
struct always_false {
  static constexpr bool value = false;
};
TEST(MathFunctions, filter_map_empty) {
  auto args = stan::math::filter_map<always_true>([](auto&& arg) { return arg; }, std::make_tuple());
  EXPECT_EQ(std::tuple_size_v<decltype(args)>, 0);
}

TEST(MathFunctions, filter_map_all_true) {
  auto args = stan::math::filter_map<always_true>([](auto&& arg) { return arg; }, std::make_tuple(1, 2, 3));
  static_assert(std::tuple_size_v<decltype(args)> == 3, "tuple size should be 3!");
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(args)>, int>));
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, decltype(args)>, int>));
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<2, decltype(args)>, int>));
}
TEST(MathFunctions, filter_map_all_false) {
  auto args = stan::math::filter_map<always_false>([](auto&& arg) { return arg; }, std::make_tuple(1, 2, 3));
  static_assert(std::tuple_size_v<decltype(args)> == 0, "tuple size should be 0!");
}

TEST(MathFunctions, filter_map_first_true) {
  auto args = stan::math::filter_map<std::is_floating_point>([](auto&& arg) { return arg; }, std::make_tuple(1.0, 2.0, 3, 4.0, 5, 6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 4, "tuple size should be 4!");
}

TEST(MathFunctions, filter_map_first_false) {
  auto args = stan::math::filter_map<std::is_floating_point>([](auto&& arg) { return arg; }, std::make_tuple(1, 2, 3, 4.0, 5, 6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 2, "tuple size should be 2!");
}

namespace internal {
template <typename T>
struct is_any_floating_point {
  static constexpr bool value = std::is_floating_point_v<T>;
};
template <typename... Types>
struct is_any_floating_point<std::tuple<Types...>> {
  static constexpr bool value = (is_any_floating_point<Types>::value || ...);
};

}
template <typename... Types>
struct is_any_floating_point {
  static constexpr bool value = (internal::is_any_floating_point<Types>::value || ...);
};

TEST(MathFunctions, filter_map_inner_tuple) {
  auto args = stan::math::filter_map<is_any_floating_point>(
     [](auto&& arg) {
      return arg;
     }, std::make_tuple(1, 2, 3, 4.0, 5, std::make_tuple(1.0, 2, 3.0), 6.0));
  static_assert(std::tuple_size_v<decltype(args)> == 3, "tuple size should be 3!");
}



}

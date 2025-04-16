#include <tuple>
#include <type_traits>
#include <string>
#include <utility>
#include <gtest/gtest.h>
#include <stan/math/prim/functor/tuple_concat.hpp>

// Test that forwarding a single tuple returns the same tuple.
TEST(TupleConcat, SingleTuple) {
  auto t = std::make_tuple(1, 2.0, 'a');
  auto result = stan::math::tuple_concat(t);
  static_assert(std::tuple_size<decltype(result)>::value == 3,
                "Result should have 3 elements");
  EXPECT_EQ(std::get<0>(result), 1);
  EXPECT_DOUBLE_EQ(std::get<1>(result), 2.0);
  EXPECT_EQ(std::get<2>(result), 'a');
}

// Test concatenating two tuples.
TEST(TupleConcat, TwoTuples) {
  auto t1 = std::make_tuple(1, 2.0);
  auto t2 = std::make_tuple('a', std::string("hello"));
  auto result = stan::math::tuple_concat(t1, t2);
  static_assert(std::tuple_size<decltype(result)>::value == 4,
                "Resulting tuple size should be 4");
  EXPECT_EQ(std::get<0>(result), 1);
  EXPECT_DOUBLE_EQ(std::get<1>(result), 2.0);
  EXPECT_EQ(std::get<2>(result), 'a');
  EXPECT_EQ(std::get<3>(result), std::string("hello"));
}

// Test concatenating three tuples.
TEST(TupleConcat, ThreeTuples) {
  auto t1 = std::make_tuple(1);
  auto t2 = std::make_tuple(2.0);
  auto t3 = std::make_tuple('a', std::string("world"));
  auto result = stan::math::tuple_concat(t1, t2, t3);
  static_assert(std::tuple_size<decltype(result)>::value == 4,
                "Resulting tuple size should be 4");
  EXPECT_EQ(std::get<0>(result), 1);
  EXPECT_DOUBLE_EQ(std::get<1>(result), 2.0);
  EXPECT_EQ(std::get<2>(result), 'a');
  EXPECT_EQ(std::get<3>(result), std::string("world"));
}

// Test concatenating more than three tuples.
TEST(TupleConcat, MultipleTuples) {
  auto t1 = std::make_tuple(1);
  auto t2 = std::make_tuple(2.0);
  auto t3 = std::make_tuple('a');
  auto t4 = std::make_tuple(std::string("test"));
  auto result = stan::math::tuple_concat(t1, t2, t3, t4);
  static_assert(std::tuple_size<decltype(result)>::value == 4,
                "Resulting tuple size should be 4");
  EXPECT_EQ(std::get<0>(result), 1);
  EXPECT_DOUBLE_EQ(std::get<1>(result), 2.0);
  EXPECT_EQ(std::get<2>(result), 'a');
  EXPECT_EQ(std::get<3>(result), std::string("test"));
}

// Test concatenation when some of the tuples are empty.
TEST(TupleConcat, EmptyTuples) {
  auto empty = std::make_tuple();
  auto t = std::make_tuple(42);
  auto result = stan::math::tuple_concat(empty, t, empty);
  static_assert(std::tuple_size<decltype(result)>::value == 1,
                "Resulting tuple size should be 1");
  EXPECT_EQ(std::get<0>(result), 42);
}

// Test that rvalue and lvalue forwarding is preserved.
TEST(TupleConcat, RvalueAndLvalueForwarding) {
  int a = 42;
  // Using forward_as_tuple to create a tuple with an lvalue reference.
  auto t1 = std::forward_as_tuple(a);
  auto t2 = std::make_tuple(3.14);
  auto result = stan::math::tuple_concat(t1, t2);
  // The first element should be an lvalue reference.
  static_assert(std::is_same_v<std::tuple_element_t<0, decltype(result)>, int&>,
                "First element should be int&");
  static_assert(
      std::is_same_v<std::tuple_element_t<1, decltype(result)>, double&>,
      "Second element should be double");
  a = 100;
  EXPECT_EQ(std::get<0>(result), 100);
}

// Test three tuples of differing sizes to ensure that all elements are included
// (this test may expose issues if the wrong tuple size is used for one of the
// inputs).
TEST(TupleConcat, ThreeTuplesDifferentSizes) {
  auto t1 = std::make_tuple(1, 2, 3);
  auto t2 = std::make_tuple(4);
  auto t3 = std::make_tuple(5, 6);  // t3 has 2 elements.
  auto result = stan::math::tuple_concat(t1, t2, t3);
  static_assert(std::tuple_size<decltype(result)>::value == 6,
                "Resulting tuple size should be 6");
  EXPECT_EQ(std::get<0>(result), 1);
  EXPECT_EQ(std::get<1>(result), 2);
  EXPECT_EQ(std::get<2>(result), 3);
  EXPECT_EQ(std::get<3>(result), 4);
  EXPECT_EQ(std::get<4>(result), 5);
  EXPECT_EQ(std::get<5>(result), 6);
}

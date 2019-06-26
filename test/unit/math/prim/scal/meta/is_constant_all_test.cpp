#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

template <typename... Ts>
void expect_is_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_TRUE(temp);
}
TEST(MetaTraits, isConstantStruct) {
  expect_is_const<>();
  expect_is_const<int>();
  expect_is_const<double>();
  expect_is_const<float>();
  expect_is_const<unsigned int>();
  expect_is_const<int32_t>();
  expect_is_const<int, int>();
  expect_is_const<double, double, double>();
  expect_is_const<float, float, float, float>();
  expect_is_const<int32_t, int32_t, int32_t, int32_t>();
}

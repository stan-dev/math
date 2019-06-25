#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

template <typename... Ts>
void expect_not_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_FALSE(temp);
}

TEST(MetaTraits, isConstantStruct) {
  expect_not_const<stan::math::var>();
  expect_not_const<stan::math::var, double>();
  expect_not_const<stan::math::var, double, stan::math::var>();
  expect_not_const<stan::math::var, double, double>();
  expect_not_const<stan::math::var, stan::math::var, stan::math::var,
                   stan::math::var, stan::math::var, double, double>();
}

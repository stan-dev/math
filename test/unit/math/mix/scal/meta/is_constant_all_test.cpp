#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

template <typename... Ts>
void expect_not_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_FALSE(temp);
}
TEST(MetaTraits, isConstant) {
  using stan::math::fvar;
  using stan::math::var;

  expect_not_const<fvar<var> >();
}

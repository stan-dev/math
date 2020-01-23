#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

template <typename... Ts>
void expect_not_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_FALSE(temp);
}

TEST(MathMetaFwd, isConstant) {
  using stan::math::fvar;

  expect_not_const<fvar<double> >();
}

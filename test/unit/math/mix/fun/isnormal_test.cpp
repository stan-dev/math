#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_isnormal() {
  using stan::math::isnormal;
  using std::isnormal;
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isnormal(inf));
  EXPECT_FALSE(isnormal(-inf));
  EXPECT_FALSE(isnormal(nan));
  EXPECT_TRUE(isnormal(T(1)));
  EXPECT_TRUE(isnormal(T(1.0)));
  EXPECT_FALSE(isnormal(T(0)));
  EXPECT_FALSE(isnormal(T(0.0)));
  EXPECT_TRUE(isnormal(T(-1)));
  EXPECT_TRUE(isnormal(T(-1.0)));
}

TEST(mixFun, isnormal) {
  using stan::math::fvar;
  using stan::math::var;
  expect_isnormal<double>();
  expect_isnormal<var>();
  expect_isnormal<fvar<double>>();
  expect_isnormal<fvar<fvar<double>>>();
  expect_isnormal<fvar<var>>();
  expect_isnormal<fvar<fvar<var>>>();
}

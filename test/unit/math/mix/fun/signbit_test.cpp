#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_signbit() {
  using stan::math::signbit;
  using std::numeric_limits;
  using std::signbit;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(signbit(inf));
  EXPECT_TRUE(signbit(-inf));
  EXPECT_FALSE(signbit(nan));
  EXPECT_FALSE(signbit(T(1)));
  EXPECT_FALSE(signbit(T(1.0)));
  EXPECT_FALSE(signbit(T(0)));
  EXPECT_FALSE(signbit(T(0.0)));
  EXPECT_TRUE(signbit(T(-0.0)));
  EXPECT_TRUE(signbit(T(-1)));
  EXPECT_TRUE(signbit(T(-1.0)));
}

TEST(mixFun, signbit) {
  using stan::math::fvar;
  using stan::math::var;
  expect_signbit<double>();
  expect_signbit<var>();
  expect_signbit<fvar<double>>();
  expect_signbit<fvar<fvar<double>>>();
  expect_signbit<fvar<var>>();
  expect_signbit<fvar<fvar<var>>>();
}

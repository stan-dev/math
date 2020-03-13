#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_isfinite() {
  using stan::math::isfinite;
  using std::isfinite;
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isfinite(inf));
  EXPECT_FALSE(isfinite(-inf));
  EXPECT_FALSE(isfinite(nan));
  EXPECT_TRUE(isfinite(T(1)));
  EXPECT_TRUE(isfinite(T(1.0)));
  EXPECT_TRUE(isfinite(T(0)));
  EXPECT_TRUE(isfinite(T(0.0)));
  EXPECT_TRUE(isfinite(T(-1)));
  EXPECT_TRUE(isfinite(T(-1.0)));
}

TEST(mixFun, isfinite) {
  using stan::math::fvar;
  using stan::math::var;
  expect_isfinite<double>();
  expect_isfinite<var>();
  expect_isfinite<fvar<double>>();
  expect_isfinite<fvar<fvar<double>>>();
  expect_isfinite<fvar<var>>();
  expect_isfinite<fvar<fvar<var>>>();
}

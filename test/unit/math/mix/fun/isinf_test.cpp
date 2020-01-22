#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_isinf() {
  using stan::math::isinf;
  using std::isinf;
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  EXPECT_TRUE(isinf(inf));
  EXPECT_TRUE(isinf(-inf));
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isinf(nan));
  EXPECT_FALSE(isinf(T(1)));
  EXPECT_FALSE(isinf(T(1.0)));
  EXPECT_FALSE(isinf(T(0)));
  EXPECT_FALSE(isinf(T(0.0)));
  EXPECT_FALSE(isinf(T(-1)));
  EXPECT_FALSE(isinf(T(-1.0)));
}

TEST(mixFun, isinf) {
  using stan::math::fvar;
  using stan::math::var;
  expect_isinf<double>();
  expect_isinf<var>();
  expect_isinf<fvar<double>>();
  expect_isinf<fvar<fvar<double>>>();
  expect_isinf<fvar<var>>();
  expect_isinf<fvar<fvar<var>>>();
}

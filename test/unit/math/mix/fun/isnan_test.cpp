#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_isnan_include_std() {
  // C++ idiom for clients of the math library
  using std::isnan;  // autodiff brought in by ADL
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isnan(inf));
  EXPECT_FALSE(isnan(-inf));
  EXPECT_TRUE(isnan(nan));
  EXPECT_FALSE(isnan(T(1)));
  EXPECT_FALSE(isnan(T(1.0)));
  EXPECT_FALSE(isnan(T(0)));
  EXPECT_FALSE(isnan(T(0.0)));
  EXPECT_FALSE(isnan(T(-1)));
  EXPECT_FALSE(isnan(T(-1.0)));
}

template <typename T>
void expect_isnan_include_stan() {
  // generated model code inclusion pattern
  using namespace stan::math;  // NOLINT(build/namespaces)
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isnan(inf));
  EXPECT_FALSE(isnan(-inf));
  EXPECT_TRUE(isnan(nan));
  EXPECT_FALSE(isnan(T(1)));
  EXPECT_FALSE(isnan(T(1.0)));
  EXPECT_FALSE(isnan(T(0)));
  EXPECT_FALSE(isnan(T(0.0)));
  EXPECT_FALSE(isnan(T(-1)));
  EXPECT_FALSE(isnan(T(-1.0)));
}

template <typename T>
void expect_isnan() {
  expect_isnan_include_std<T>();
  expect_isnan_include_stan<T>();
}

TEST(mixFun, isnan) {
  using stan::math::fvar;
  using stan::math::var;
  expect_isnan<double>();
  expect_isnan<var>();
  expect_isnan<fvar<double>>();
  expect_isnan<fvar<fvar<double>>>();
  expect_isnan<fvar<var>>();
  expect_isnan<fvar<fvar<var>>>();
}

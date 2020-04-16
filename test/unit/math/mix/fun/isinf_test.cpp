#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename T>
void expect_isinf_include_std() {
  // standard usage:
  // - std::isinf explicitly brought in,
  // - stan::math::isinf brought in by argument-dependent lookup (ADL)
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

template <typename T>
void expect_isinf_include_stan() {
  // usage in generated model code:
  // - stan::math::isinf explicitly brought in
  using namespace stan::math;
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

template <typename T>
void expect_isinf() {
  expect_isinf_include_std<T>();
  expect_isinf_include_stan<T>();
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

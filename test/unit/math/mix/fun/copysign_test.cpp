#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

template <typename T>
void expect_copysign() {
  using stan::math::copysign;
  using stan::math::signbit;
  using std::copysign;
  using std::numeric_limits;
  using std::signbit;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<T> ys{inf, -inf, -1, 0, 1};  // inf, -inf, -1.0, 0.0, 1.0 };

  for (const T& x : ys) {
    for (const T& y : ys) {
      EXPECT_EQ(signbit(y), signbit(copysign(T(x), y)));
      EXPECT_EQ(signbit(y), signbit(copysign(x, T(y))));
      EXPECT_EQ(signbit(y), signbit(copysign(T(x), T(y))));
    }
  }
}

TEST(mixFun, copysign) {
  using stan::math::fvar;
  using stan::math::var;
  expect_copysign<double>();
  expect_copysign<var>();
  expect_copysign<fvar<double>>();
  expect_copysign<fvar<fvar<double>>>();
  expect_copysign<fvar<var>>();
  expect_copysign<fvar<fvar<var>>>();
}

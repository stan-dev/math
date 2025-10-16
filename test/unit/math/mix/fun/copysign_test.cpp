#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

template <typename U, typename V>
void expect_eq_signbit(const U& u, const V& v) {
  using stan::math::signbit;
  using std::signbit;
  EXPECT_EQ(signbit(u), signbit(v));
}

template <typename T>
void expect_copysign() {
  using stan::math::copysign;
  using stan::math::copysign_non_zero;
  using stan::math::is_nan;
  using stan::math::value_of_rec;
  using std::copysign;
  using std::numeric_limits;

  double nan = numeric_limits<double>::quiet_NaN();
  double inf = numeric_limits<double>::infinity();
  std::vector<double> vs{-inf, -1, -0.0, 0.0, 1, inf, nan};

  // real
  for (double x : vs) {
    for (double y : vs) {
      expect_eq_signbit(y, copysign(T(x), y));
      expect_eq_signbit(y, copysign(x, T(y)));
      expect_eq_signbit(y, copysign(T(x), T(y)));
    }
  }

  // complex
  for (double re1 : vs) {
    for (double im1 : vs) {
      auto x_d = std::complex<double>(re1, im1);
      auto x_t = std::complex<T>(re1, im1);
      for (double re2 : vs) {
        for (double im2 : vs) {
          auto y_d = std::complex<double>(re2, im2);
          auto y_t = std::complex<T>(re2, im2);
          if (re1 == 0) {
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_d, y_d).real()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_d, y_t).real()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_t, y_d).real()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_t, y_t).real()));
          } else {
            expect_eq_signbit(re2, copysign(x_d, y_d).real());
            expect_eq_signbit(re2, copysign(x_d, y_t).real());
            expect_eq_signbit(re2, copysign(x_t, y_d).real());
            expect_eq_signbit(re2, copysign(x_t, y_t).real());
          }
          if (im1 == 0) {
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_d, y_d).imag()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_d, y_t).imag()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_t, y_d).imag()));
            EXPECT_FLOAT_EQ(0.0, value_of_rec(copysign(x_t, y_t).imag()));
          } else {
            expect_eq_signbit(im2, copysign(x_d, y_d).imag());
            expect_eq_signbit(im2, copysign(x_d, y_t).imag());
            expect_eq_signbit(im2, copysign(x_t, y_d).imag());
            expect_eq_signbit(im2, copysign(x_t, y_t).imag());
          }
        }
      }
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

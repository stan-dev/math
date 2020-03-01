#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <complex>
#include <limits>
#include <vector>

template <typename U, typename V>
void expect_eq_signbit(const U& u, const V& v) {
  EXPECT_EQ(signbit(u), signbit(v));
}

template <typename T>
void expect_copysign() {
  using stan::math::copysign;
  using stan::math::signbit;
  using std::copysign;
  using std::numeric_limits;
  using std::signbit;

  double inf = numeric_limits<double>::infinity();
  std::vector<double> ys{inf, -inf, -1, 0, 1};

  // real
  for (double x : ys) {
    for (double y : ys) {
      expect_eq_signbit(y, copysign(T(x), y));
      expect_eq_signbit(y, copysign(x, T(y)));
      expect_eq_signbit(y, copysign(T(x), T(y)));
    }
  }
  // complex
  for (double re1 : ys) {
    for (double im1 : ys) {
      auto x_d = std::complex<double>(re1, im1);
      auto x_t = std::complex<T>(re1, im1);
      for (double re2 : ys) {
        for (double im2 : ys) {
          auto y_d = std::complex<double>(re2, im2);
          auto y_t = std::complex<T>(re2, im2);
          expect_eq_signbit(re2, copysign(x_d, y_d).real());
          expect_eq_signbit(re2, copysign(x_d, y_t).real());
          expect_eq_signbit(re2, copysign(x_t, y_d).real());
          expect_eq_signbit(re2, copysign(x_t, y_t).real());
          expect_eq_signbit(im2, copysign(x_d, y_d).imag());
          expect_eq_signbit(im2, copysign(x_d, y_t).imag());
          expect_eq_signbit(im2, copysign(x_t, y_d).imag());
          expect_eq_signbit(im2, copysign(x_t, y_t).imag());
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

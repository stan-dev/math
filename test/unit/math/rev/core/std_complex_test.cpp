#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>

// naughty but convenient for test scope
typedef stan::math::var var_t;
typedef std::complex<stan::math::var> cvar_t;
typedef std::complex<double> cdouble_t;

void expect_complex(double re, double im, const cvar_t& y) {
  EXPECT_DOUBLE_EQ(re, y.real().val());
  EXPECT_DOUBLE_EQ(im, y.imag().val());
}
void expect_complex(const cdouble_t& x, const cvar_t& y) {
  expect_complex(x.real(), x.imag(), y);
}

TEST(mathRevCore, stdComplexConstructor1) {
  // constructor (1), no defaults
  var_t x = 1;
  var_t y = 2;
  cvar_t z{x, y};
  expect_complex(1, 2, z);

  // constructor (1), default imaginary
  cvar_t a{x};
  expect_complex(1, 0, a);

  // constructor (1), default real and imaginary
  cvar_t b{};
  expect_complex(0, 0, b);

  // direct double and mixed
  cvar_t c{x, 2.0};
  expect_complex(1, 2, c);

  cvar_t d{1.0, y};
  expect_complex(1, 2, d);

  cvar_t e{1.0, 2.0};
  expect_complex(1, 2, e);

  // int and mixed
  cvar_t f{x, 2};
  expect_complex(1, 2, f);

  cvar_t g{1, y};
  expect_complex(1, 2, g);

  cvar_t h{1, 2};
  expect_complex(1, 2, h);

  cvar_t j{1.0, 2};
  expect_complex(1, 2, j);

  cvar_t k{1, 2.0};
  expect_complex(1, 2, k);

  cvar_t l{1.0, 2.0};
  expect_complex(1, 2, l);
}
TEST(mathRevCore, stdComplexConstructor2) {
  var_t x = 1;
  var_t y = 2;
  cvar_t a{x, y};
  cvar_t b{a};
  // require pimpl equality
  EXPECT_EQ(a.real().vi_, b.real().vi_);
  EXPECT_EQ(a.imag().vi_, b.imag().vi_);
}
TEST(mathRevCore, stdComplexConstructor3) {
  cdouble_t a(1, 2);
  cvar_t b(a);
  expect_complex(1, 2, b);
}
TEST(MathRevCore, stdComplexReal1) {
  cvar_t a(3, -1);
  EXPECT_DOUBLE_EQ(3, a.real().val());
}
TEST(MathRevCore, stdComplexReal2) {
  cvar_t a(3, -1);
  a.real(2.7);
  EXPECT_DOUBLE_EQ(2.7, a.real().val());
  EXPECT_DOUBLE_EQ(-1, a.imag().val());
}
TEST(MathRevCore, stdComplexImag1) {
  cvar_t a(3, -1);
  EXPECT_DOUBLE_EQ(-1, a.imag().val());
}
TEST(MathRevCore, stdComplexImag2) {
  cvar_t a(3, -1);
  a.imag(2.7);
  EXPECT_DOUBLE_EQ(2.7, a.imag().val());
  EXPECT_DOUBLE_EQ(3, a.real().val());
}
TEST(MathRevCore, stdComplexOperatorEquals1) {
  cvar_t a(1, 2);
  var_t b = 3;
  a = b;
  // require pimpl equality
  EXPECT_EQ(b.vi_, a.real().vi_);
  // require imaginary value to be 0
  EXPECT_EQ(0, a.imag().val());
  // require return of *this
  var_t c = 4;
  auto ptr1 = &a;
  auto ptr2 = &(a = c);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorEquals2) {
  cvar_t a(1, 2);
  cvar_t b(3, 4);
  a = b;
  // require pimpl equality
  EXPECT_EQ(b.real().vi_, a.real().vi_);
  EXPECT_EQ(b.imag().vi_, a.imag().vi_);
  // require return of *this
  cvar_t c(5, 6);
  auto ptr1 = &a;
  auto ptr2 = &(a = c);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorEquals3) {
  cvar_t a(1, 2);
  cdouble_t b(3, 4);
  a = b;
  EXPECT_EQ(3, a.real().val());
  EXPECT_EQ(4, a.imag().val());
  // require return of *this
  cdouble_t c(5, 6);
  auto ptr1 = &a;
  auto ptr2 = &(a = c);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorAddEquals1) {
  cdouble_t ad(1, 2);
  ad += 3;

  cvar_t a(1, 2);
  var_t b = 3;
  a += b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a += b);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorSubtractEquals2) {
  cdouble_t ad(1, 2);
  ad -= 3;

  cvar_t a(1, 2);
  var_t b = 3;
  a -= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a -= b);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorMultiplyEquals3) {
  cdouble_t ad(2, 5);
  ad *= 3;

  cvar_t a(2, 5);
  var_t b = 3;
  a *= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a *= b);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorDivideEquals4) {
  cdouble_t ad(2, 5);
  ad /= 3;

  cvar_t a(2, 5);
  var_t b = 3;
  a /= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a /= b);
  EXPECT_EQ(ptr1, ptr2);
}
TEST(MathRevCore, stdComplexOperatorAddEquals5) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 4);
  ad += bd;

  cvar_t a(1, 2);
  cvar_t b(3, 4);
  a += b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a += b);
  EXPECT_EQ(ptr1, ptr2);

  cvar_t c(1, 2);
  cdouble_t d(3, 4);
  c += d;
  expect_complex(ad, c);
}
TEST(MathRevCore, stdComplexOperatorSubtractEquals6) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 4);
  ad -= bd;

  cvar_t a(1, 2);
  cvar_t b(3, 4);
  a -= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a += b);
  EXPECT_EQ(ptr1, ptr2);

  cvar_t c(1, 2);
  cdouble_t d(3, 4);
  c -= d;
  expect_complex(ad, c);
}
TEST(MathRevCore, stdComplexOperatorMultiplyEquals7) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 4);
  ad *= bd;

  cvar_t a(1, 2);
  cvar_t b(3, 4);
  a *= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a *= b);
  EXPECT_EQ(ptr1, ptr2);

  cvar_t c(1, 2);
  cdouble_t d(3, 4);
  c *= d;
  expect_complex(ad, c);
}
TEST(MathRevCore, stdComplexOperatorDivideEquals8) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 4);
  ad /= bd;

  cvar_t a(1, 2);
  cvar_t b(3, 4);
  a /= b;
  expect_complex(ad, a);

  auto ptr1 = &a;
  auto ptr2 = &(a /= b);
  EXPECT_EQ(ptr1, ptr2);

  cvar_t c(1, 2);
  cdouble_t d(3, 4);
  c /= d;
  expect_complex(ad, c);
}
TEST(MathRevCore, stdComplexOperatorUnaryPlus1) {
  cdouble_t ad(1, 2);
  cdouble_t bd = +ad;
  cvar_t a(1, 2);
  cvar_t b = +a;
  expect_complex(bd, b);
  // expect pointer equality before and after
  EXPECT_EQ(a.real().vi_, b.real().vi_);
}
TEST(MathRevCore, stdComplexOperatorUnaryNegative2) {
  cdouble_t ad(1, 2);
  cdouble_t bd = -ad;
  cvar_t a(1, 2);
  cvar_t b = -a;
  expect_complex(bd, b);
}

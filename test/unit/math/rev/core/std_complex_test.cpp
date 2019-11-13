#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <sstream>
#include <vector>

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
TEST(mathRevCore, stdComplexReal1) {
  cvar_t a(3, -1);
  EXPECT_DOUBLE_EQ(3, a.real().val());
}
TEST(mathRevCore, stdComplexReal2) {
  cvar_t a(3, -1);
  a.real(2.7);
  EXPECT_DOUBLE_EQ(2.7, a.real().val());
  EXPECT_DOUBLE_EQ(-1, a.imag().val());
}
TEST(mathRevCore, stdComplexImag1) {
  cvar_t a(3, -1);
  EXPECT_DOUBLE_EQ(-1, a.imag().val());
}
TEST(mathRevCore, stdComplexImag2) {
  cvar_t a(3, -1);
  a.imag(2.7);
  EXPECT_DOUBLE_EQ(2.7, a.imag().val());
  EXPECT_DOUBLE_EQ(3, a.real().val());
}
TEST(mathRevCore, stdComplexOperatorEquals1) {
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
TEST(mathRevCore, stdComplexOperatorEquals2) {
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
TEST(mathRevCore, stdComplexOperatorEquals3) {
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
TEST(mathRevCore, stdComplexOperatorAddEquals1) {
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
TEST(mathRevCore, stdComplexOperatorSubtractEquals2) {
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
TEST(mathRevCore, stdComplexOperatorMultiplyEquals3) {
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
TEST(mathRevCore, stdComplexOperatorDivideEquals4) {
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
TEST(mathRevCore, stdComplexOperatorAddEquals5) {
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
TEST(mathRevCore, stdComplexOperatorSubtractEquals6) {
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
TEST(mathRevCore, stdComplexOperatorMultiplyEquals7) {
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
TEST(mathRevCore, stdComplexOperatorDivideEquals8) {
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
TEST(mathRevCore, stdComplexOperatorUnaryPlus1) {
  cdouble_t ad(1, 2);
  cdouble_t bd = +ad;
  cvar_t a(1, 2);
  cvar_t b = +a;
  expect_complex(bd, b);
  // expect pointer equality before and after
  EXPECT_EQ(a.real().vi_, b.real().vi_);
}
TEST(mathRevCore, stdComplexOperatorUnaryNegative2) {
  cdouble_t ad(1, 2);
  cdouble_t bd = -ad;
  cvar_t a(1, 2);
  cvar_t b = -a;
  expect_complex(bd, b);
}
TEST(mathRevCore, stdComplexOperatorAdd1) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  cdouble_t cd = ad + bd;
  cvar_t a(1, 2);
  cvar_t b(3, 7);
  cvar_t c = a + b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorAdd2) {
  cdouble_t ad(1, 2);
  double bd = 3;
  cdouble_t cd = ad + bd;
  cvar_t a(1, 2);
  var_t b = 3;
  cvar_t c = a + b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorAdd3) {
  double ad = 1;
  cdouble_t bd(3, 7);
  cdouble_t cd = ad + bd;
  var_t a = 1;
  cvar_t b(3, 7);
  cvar_t c = a + b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorSubtract4) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  cdouble_t cd = ad - bd;
  cvar_t a(1, 2);
  cvar_t b(3, 7);
  cvar_t c = a - b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorSubtract5) {
  cdouble_t ad(1, 2);
  double bd = 3;
  cdouble_t cd = ad - bd;
  cvar_t a(1, 2);
  var_t b = 3;
  cvar_t c = a - b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorSubtract6) {
  double ad = 1;
  cdouble_t bd(3, 7);
  cdouble_t cd = ad - bd;
  var_t a = 1;
  cvar_t b(3, 7);
  cvar_t c = a - b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorMultiply7) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  cdouble_t cd = ad * bd;
  cvar_t a(1, 2);
  cvar_t b(3, 7);
  cvar_t c = a * b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorMultiply8) {
  cdouble_t ad(1, 2);
  double bd = 3;
  cdouble_t cd = ad * bd;
  cvar_t a(1, 2);
  var_t b = 3;
  cvar_t c = a * b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorMultiply9) {
  double ad = 1;
  cdouble_t bd(3, 7);
  cdouble_t cd = ad * bd;
  var_t a = 1;
  cvar_t b(3, 7);
  cvar_t c = a * b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorDivide10) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  cdouble_t cd = ad / bd;
  cvar_t a(1, 2);
  cvar_t b(3, 7);
  cvar_t c = a / b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorDivide11) {
  cdouble_t ad(1, 2);
  double bd = 3;
  cdouble_t cd = ad / bd;
  cvar_t a(1, 2);
  var_t b = 3;
  cvar_t c = a / b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorDivide12) {
  double ad = 1;
  cdouble_t bd(3, 7);
  cdouble_t cd = ad / bd;
  var_t a = 1;
  cvar_t b(3, 7);
  cvar_t c = a / b;
  expect_complex(cd, c);
}
TEST(mathRevCore, stdComplexOperatorCompare1) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  bool cd = ad == bd;
  bool dd = ad == ad;

  cvar_t a(1, 2);
  cvar_t b(3, 7);
  bool c = a == b;
  bool d = a == a;
  EXPECT_EQ(cd, c);
  EXPECT_EQ(dd, d);
}
TEST(mathRevCore, stdComplexOperatorCompare2) {
  cdouble_t ad(1, 0);
  cdouble_t bd(2, 3);
  double cd = 1;
  double dd = 3.2;

  bool ed = ad == cd;
  bool fd = ad == dd;
  bool gd = bd == cd;
  bool hd = bd == dd;

  cvar_t a(1, 0);
  cvar_t b(2, 3);
  var_t c = 1;
  var_t d = 3.2;

  bool e = a == c;
  bool f = a == d;
  bool g = b == c;
  bool h = b == d;
  EXPECT_EQ(ed, e);
  EXPECT_EQ(fd, f);
  EXPECT_EQ(gd, g);
  EXPECT_EQ(hd, h);
}
TEST(mathRevCore, stdComplexOperatorCompare3) {
  cdouble_t ad(1, 0);
  cdouble_t bd(2, 3);
  double cd = 1;
  double dd = 3.2;

  bool ed = cd == ad;
  bool fd = cd == bd;
  bool gd = dd == ad;
  bool hd = dd == bd;

  cvar_t a(1, 0);
  cvar_t b(2, 3);
  var_t c = 1;
  var_t d = 3.2;

  bool e = c == a;
  bool f = c == b;
  bool g = d == a;
  bool h = d == b;
  EXPECT_EQ(ed, e);
  EXPECT_EQ(fd, f);
  EXPECT_EQ(gd, g);
  EXPECT_EQ(hd, h);
}
TEST(mathRevCore, stdComplexOperatorCompareUneq4) {
  cdouble_t ad(1, 2);
  cdouble_t bd(3, 7);
  bool cd = ad != bd;
  bool dd = ad != ad;

  cvar_t a(1, 2);
  cvar_t b(3, 7);
  bool c = a != b;
  bool d = a != a;
  EXPECT_EQ(cd, c);
  EXPECT_EQ(dd, d);
}
TEST(mathRevCore, stdComplexOperatorCompareUneq5) {
  cdouble_t ad(1, 0);
  cdouble_t bd(2, 3);
  double cd = 1;
  double dd = 3.2;

  bool ed = ad != cd;
  bool fd = ad != dd;
  bool gd = bd != cd;
  bool hd = bd != dd;

  cvar_t a(1, 0);
  cvar_t b(2, 3);
  var_t c = 1;
  var_t d = 3.2;

  bool e = a != c;
  bool f = a != d;
  bool g = b != c;
  bool h = b != d;
  EXPECT_EQ(ed, e);
  EXPECT_EQ(fd, f);
  EXPECT_EQ(gd, g);
  EXPECT_EQ(hd, h);
}
TEST(mathRevCore, stdComplexOperatorCompareUneq6) {
  cdouble_t ad(1, 0);
  cdouble_t bd(2, 3);
  double cd = 1;
  double dd = 3.2;

  bool ed = cd != ad;
  bool fd = cd != bd;
  bool gd = dd != ad;
  bool hd = dd != bd;

  cvar_t a(1, 0);
  cvar_t b(2, 3);
  var_t c = 1;
  var_t d = 3.2;

  bool e = c != a;
  bool f = c != b;
  bool g = d != a;
  bool h = d != b;
  EXPECT_EQ(ed, e);
  EXPECT_EQ(fd, f);
  EXPECT_EQ(gd, g);
  EXPECT_EQ(hd, h);
}
// the operator streams use the compiler-supplied implementations
// this is generally unspecified behavior, but seems to work
TEST(mathRevCore, stanMathOperatorStreamOut1) {
  cdouble_t ad(1, 2);
  cvar_t a(1, 2);

  std::stringstream ssd;
  ssd << ad;
  std::string sd = ssd.str();
  std::stringstream ss;
  ss << a;
  std::string s = ss.str();
  EXPECT_EQ(sd, s);
}
TEST(mathRevCore, stanMathOperatorStreamIn2) {
  std::stringstream s1;
  s1 << "(1, 2)";
  cvar_t a;
  s1 >> a;
  expect_complex(1, 2, a);

  std::stringstream s2;
  s2 << "(1)";
  s2 >> a;
  expect_complex(1, 0, a);

  std::stringstream s3;
  s3 << "1";
  s3 >> a;
  expect_complex(1, 0, a);
}
TEST(mathRevCore, stdRealExternal1) {
  var_t a = 1;
  var_t b = 2;
  cvar_t c(a, b);
  var_t ca = std::real(c);
  EXPECT_EQ(a.vi_, ca.vi_);
}
TEST(mathRevCore, stdImagExternal1) {
  var_t a = 1;
  var_t b = 2;
  cvar_t c(a, b);
  var_t cb = std::imag(c);
  EXPECT_EQ(b.vi_, cb.vi_);
}
TEST(mathRevCore, stdAbsExternal1) {
  cdouble_t ad(1, 2);
  double bd = std::abs(ad);
  cvar_t a(1, 2);
  var_t b = std::abs(a);
  EXPECT_DOUBLE_EQ(bd, b.val());
}
TEST(mathRevCore, stdArgExternal1) {
  cdouble_t ad(1, 2);
  double bd = std::arg(ad);
  cvar_t a(1, 2);
  var_t b = std::arg(a);
  EXPECT_DOUBLE_EQ(bd, b.val());
}
TEST(mathRevCore, stdNormExternal1) {
  cdouble_t ad(1, 2);
  double bd = std::norm(ad);
  cvar_t a(1, 2);
  var_t b = std::norm(a);
  EXPECT_DOUBLE_EQ(bd, b.val());
}
TEST(mathRevCore, stdSquareConj1) {
  cdouble_t ad(1, 2);
  cdouble_t bd = std::conj(ad);
  cvar_t a(1, 2);
  cvar_t b = std::conj(a);
  expect_complex(bd, b);
}
TEST(mathRevCore, stdSquareProj1) {
  double inf = std::numeric_limits<double>::infinity();
  std::vector<double> args{-1, 0, 1, inf, -inf};
  for (double re : args) {
    for (double im : args) {
      cdouble_t ad(re, im);
      cdouble_t bd = std::proj(ad);
      cvar_t a(re, im);
      cvar_t b = std::proj(a);
      expect_complex(bd, b);
    }
  }
}
TEST(mathRevCore, stdPolar1) {
  double r_d = 0.5;
  double theta_d = 1.3;
  cdouble_t a_d = std::polar(r_d, theta_d);

  var_t r = 0.5;
  var_t theta = 1.3;
  cvar_t a = std::polar(r, theta);
  expect_complex(a_d, a);

  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  // these are undefined behaviors, but we return NaN
  var_t neg1_v = -1;
  var_t one_v = 1;
  var_t two_v = 2;
  var_t nan_v = nan;
  var_t inf_v = inf;
  EXPECT_TRUE(stan::math::is_nan(std::polar(neg1_v, two_v).real()));
  EXPECT_TRUE(stan::math::is_nan(std::polar(nan_v, two_v).real()));
  EXPECT_TRUE(stan::math::is_nan(std::polar(one_v, inf_v).real()));
}

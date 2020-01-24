#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// DAUNTING COMBINATORICS FOR BINARY OPERATIONS
// (7 real)     i,  d,  v,  fd,  ffd,  fv,  ffv
// (6 complex)  -,  cd, cv, cfd, cffd, cfv, cffv
// ----------
// 13 * 13 = 169 test combinations!  uh oh

using mvar_t = Eigen::Matrix<stan::math::var, -1, -1>;

using var_t = stan::math::var;
using fvar_d_t = stan::math::fvar<double>;
using fvar_fvar_d_t = stan::math::fvar<fvar_d_t>;
using fvar_v_t = stan::math::fvar<stan::math::var>;
using fvar_fvar_v_t = stan::math::fvar<fvar_v_t>;

using cdouble_t = std::complex<double>;
using cvar_t = std::complex<stan::math::var>;
using cfvar_d_t = std::complex<fvar_d_t>;
using cfvar_fvar_d_t = std::complex<fvar_fvar_d_t>;
using cfvar_v_t = std::complex<fvar_v_t>;
using cfvar_fvar_v_t = std::complex<fvar_fvar_v_t>;

template <typename T>
void expect_identity_matrix(const T& I) {
  for (int i = 0; i < I.rows(); ++i) {
    EXPECT_NEAR(1, value_of_rec(I(i, i)), 1e-8);
    for (int j = 0; j < i; ++j) {
      EXPECT_NEAR(0, value_of_rec(I(i, j)), 1e-8);
      EXPECT_NEAR(0, value_of_rec(I(j, i)), 1e-8);
    }
  }
}

template <typename F>
void expect_reduction(const F& f) {
  cvar_t x(1, 2);
  var_t y = f(x);

  cdouble_t xd(1, 2);
  double yd = f(xd);
  EXPECT_FLOAT_EQ(yd, y.val());
}

void expect_complex(double re, double im, const cdouble_t& z) {
  stan::test::expect_near_rel("complex real value", re, z.real());
  stan::test::expect_near_rel("complex imag value", im, z.imag());
}
template <typename T>
void expect_complex(double re, double im, const std::complex<T>& z) {
  auto z_reduced
      = std::complex<decltype(z.real().val())>{z.real().val(), z.imag().val()};
  expect_complex(re, im, z_reduced);
}

template <typename T>
void expect_complex(const cdouble_t& x, const std::complex<T>& y) {
  expect_complex(x.real(), x.imag(), y);
}

std::vector<double> common_non_neg_vals() {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double pos_zero = 0.0;
  double neg_zero = -0.0;
  return {0.0, 1.3, 2.1};
}

std::vector<double> common_vals() {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double pos_zero = 0.0;
  double neg_zero = -0.0;
  return {-4, -2.5, -1.5, -0.3, -0.0, 0.0, 1.3, 2.1, 3.9};
}

template <typename F>
void expect_complex_common(const F& f) {
  for (double a : common_vals()) {
    for (double b : common_vals()) {
      expect_complex(f(cdouble_t{a, b}), f(cvar_t{a, b}));
    }
  }
}

template <typename F>
void expect_complex_common_binary(const F& f) {
  using stan::math::pow;
  for (double x1 : common_non_neg_vals()) {
    for (double y1 : common_non_neg_vals()) {
      for (double x2 : common_non_neg_vals()) {
        expect_complex(f(cdouble_t{x1, y1}, x2), f(cvar_t{x1, y1}, var_t{x2}));
        expect_complex(f(cdouble_t{x1, y1}, x2), f(cvar_t{x1, y1}, x2));
        expect_complex(f(x2, cdouble_t{x1, y1}), f(var_t{x2}, cvar_t{x1, y1}));
        expect_complex(f(x2, cdouble_t{x1, y1}), f(x2, cvar_t{x1, y1}));
        for (double y2 : common_non_neg_vals()) {
          expect_complex(f(cdouble_t{x1, y1}, cdouble_t{x2, y2}),
                         f(cvar_t{x1, y1}, cvar_t{x2, y2}));
        }
      }
    }
  }
}

template <typename F>
void expect_compound_assign_operator(const F& f) {
  // cvar += var
  cdouble_t ad(1, 2);
  f(ad, 3);
  cvar_t a(1, 2);
  var_t b = 3;
  f(a, b);
  expect_complex(ad, a);
  auto ptr1 = &a;
  auto ptr2 = &(a += b);
  EXPECT_EQ(ptr1, ptr2);

  // cvar += cdouble
  cvar_t e(1, 2);
  cdouble_t ed(1, 2);
  cdouble_t bd(3, 4);
  f(ed, bd);
  f(e, bd);
  expect_complex(ed, e);
  ptr1 = &e;
  ptr2 = &(e += bd);
  EXPECT_EQ(ptr1, ptr2);

  // cvar += cvar
  cdouble_t gd(1, 2);
  cdouble_t hd(3, 4);
  f(gd, hd);
  cvar_t g(1, 2);
  cvar_t h(3, 4);
  f(g, h);
  expect_complex(gd, g);
  ptr1 = &g;
  ptr2 = &(g += h);
  EXPECT_EQ(ptr1, ptr2);

  // cvar += double
  cdouble_t kd(1, 2);
  f(kd, 3.0);
  cvar_t k(1, 2);
  f(k, 3.0);
  expect_complex(kd, k);
  ptr1 = &k;
  ptr2 = &(k += 3.0);
  EXPECT_EQ(ptr1, ptr2);

  // cvar += int
  cdouble_t jd(1, 2);
  f(jd, 3);
  cvar_t j(1, 2);
  f(j, 3);
  expect_complex(jd, j);
  ptr1 = &j;
  ptr2 = &(j += 3);
  EXPECT_EQ(ptr1, ptr2);
}

template <typename T>
std::vector<T> to_array(const std::complex<T>& a) {
  return {a.real(), a.imag()};
}
template <typename T>
std::complex<T> from_array(const std::vector<T>& a) {
  return {a[0], a[1]};
}

std::vector<std::vector<double>> common_complex() {
  std::vector<std::vector<double>> zs;
  for (double x = -1; x <= 2.1; ++x)
    for (double y = -1; y <= 2.1; ++y)
      zs.push_back(std::vector<double>{x, y});
  return zs;
}

// BEGIN PR 1
// ========================================================

TEST(mathRevCore, stdIteratorTraits) {
  using itv_t = std::iterator_traits<var_t>;
  var_t v = 1.3;

  itv_t::value_type v2 = v;
  EXPECT_FLOAT_EQ(1.3, v2.val());

  itv_t::pointer v_ptr = &v;
  EXPECT_FLOAT_EQ(1.3, v_ptr->val());
  *v_ptr += 2.1;
  EXPECT_FLOAT_EQ(3.4, v.val());

  itv_t::reference v_ref = v;
  v_ref = 9.7;
  EXPECT_FLOAT_EQ(9.7, v.val());

  itv_t::difference_type a = v_ptr - v_ptr;
}

TEST(mathFwdCore, stdIteratorTraits) {
  using itv_t = std::iterator_traits<fvar_d_t>;
  fvar_d_t v = 1.3;

  itv_t::value_type v2 = v;
  EXPECT_FLOAT_EQ(1.3, v2.val());

  itv_t::pointer v_ptr = &v;
  EXPECT_FLOAT_EQ(1.3, v_ptr->val());
  *v_ptr += 2.1;
  EXPECT_FLOAT_EQ(3.4, v.val());

  itv_t::reference v_ref = v;
  v_ref = 9.7;
  EXPECT_FLOAT_EQ(9.7, v.val());

  itv_t::difference_type a = v_ptr - v_ptr;
}

template <typename T>
void expectSignBit() {
  EXPECT_TRUE(signbit(-std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(signbit(std::numeric_limits<T>::infinity()));
  EXPECT_TRUE(signbit(T{-2}));
  EXPECT_FALSE(signbit(T{0}));
  EXPECT_FALSE(signbit(T{3.9}));
}
TEST(mathFwdCore, signbit) {
  expectSignBit<var_t>();
  expectSignBit<fvar_d_t>();
  expectSignBit<fvar_fvar_d_t>();
  expectSignBit<fvar_v_t>();
  expectSignBit<fvar_fvar_v_t>();
}

template <typename T>
void expectIsInf() {
  EXPECT_FALSE(isinf(std::numeric_limits<T>::quiet_NaN()));
  EXPECT_TRUE(isinf(std::numeric_limits<T>::infinity()));
  EXPECT_TRUE(isinf(-std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isinf(T{-1}));
  EXPECT_FALSE(isinf(T{0}));
  EXPECT_FALSE(isinf(T{3e27}));
}
TEST(mathFwdCore, isinf) {
  expectIsInf<var_t>();
  expectIsInf<fvar_d_t>();
  expectIsInf<fvar_fvar_d_t>();
  expectIsInf<fvar_v_t>();
  expectIsInf<fvar_fvar_v_t>();
}

template <typename T>
void expectIsFinite() {
  EXPECT_FALSE(isfinite(std::numeric_limits<T>::quiet_NaN()));
  EXPECT_FALSE(isfinite(std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isfinite(-std::numeric_limits<T>::infinity()));
  EXPECT_TRUE(isfinite(T{-1}));
  EXPECT_TRUE(isfinite(T{0}));
  EXPECT_TRUE(isfinite(T{3e27}));
}
TEST(mathFwdCore, isfinite) {
  expectIsFinite<var_t>();
  expectIsFinite<fvar_d_t>();
  expectIsFinite<fvar_fvar_d_t>();
  expectIsFinite<fvar_v_t>();
  expectIsFinite<fvar_fvar_v_t>();
}

template <typename T>
void expectIsNan() {
  EXPECT_TRUE(isnan(std::numeric_limits<T>::quiet_NaN()));
  EXPECT_FALSE(isnan(std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isnan(-std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isnan(T{-1}));
  EXPECT_FALSE(isnan(T{0}));
  EXPECT_FALSE(isnan(T{3e27}));
}
TEST(mathFwdCore, isnan) {
  expectIsNan<var_t>();
  expectIsNan<fvar_d_t>();
  expectIsNan<fvar_fvar_d_t>();
  expectIsNan<fvar_v_t>();
  expectIsNan<fvar_fvar_v_t>();
}

template <typename T>
void expectIsNormal() {
  EXPECT_FALSE(isnormal(std::numeric_limits<T>::quiet_NaN()));
  EXPECT_FALSE(isnormal(std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isnormal(-std::numeric_limits<T>::infinity()));
  EXPECT_FALSE(isnormal(T{0}));
  EXPECT_TRUE(isnormal(T{-1}));
  EXPECT_TRUE(isnormal(T{3e27}));
}
TEST(mathFwdCore, isnormal) {
  expectIsNormal<var_t>();
  expectIsNormal<fvar_d_t>();
  expectIsNormal<fvar_fvar_d_t>();
  expectIsNormal<fvar_v_t>();
  expectIsNormal<fvar_fvar_v_t>();
}

TEST(mathFwdCore, copysignScalar) {
  auto f = [](const auto& x, const auto& y) {
    using std::copysign;
    return copysign(x, y);
  };
  std::vector<double> vals{-2.3, -1, 1.7};
  for (double x : vals) {
    for (double y : vals) {
      stan::test::expect_ad(f, x, y);
    }
  }
}

// BEGIN PR 2
// ==============================================

template <typename T>
void expect_value_of() {
  cdouble_t zd{2.0, 3.0};
  std::complex<T> z{2.0, 3.0};
  EXPECT_TRUE(zd == stan::math::value_of_rec(zd));
}

TEST(mathMix, valueOfRec) {
  expect_value_of<double>();
  expect_value_of<var_t>();
  expect_value_of<fvar_d_t>();
  expect_value_of<fvar_fvar_d_t>();
  expect_value_of<fvar_v_t>();
  expect_value_of<fvar_fvar_v_t>();
}

TEST(mathMix, copysignComplex) {
  auto f = [](const auto& x, const auto& y) {
    using stan::math::copysign;
    return to_array(copysign(from_array(x), from_array(y)));
  };
  std::vector<double> vals{-2.3, -1, 1.7};
  for (double x1 : vals) {
    for (double y1 : vals) {
      std::vector<double> z1{x1, y1};
      for (double x2 : vals) {
        for (double y2 : vals) {
          std::vector<double> z2{x2, y2};
          stan::test::expect_ad(f, z1, z2);
        }
      }
    }
  }
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

  // verifies std::complex<var>() produces (0, 0)
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
  EXPECT_FLOAT_EQ(3, a.real().val());
}
TEST(mathRevCore, stdComplexReal2) {
  // test double, int, var args
  cvar_t a(3, -1);
  a.real(2.7);
  EXPECT_FLOAT_EQ(2.7, a.real().val());
  EXPECT_FLOAT_EQ(-1, a.imag().val());

  a.real(2);
  EXPECT_FLOAT_EQ(2, a.real().val());
  EXPECT_FLOAT_EQ(-1, a.imag().val());

  var_t c(5);
  a.real(c);
  EXPECT_EQ(c.vi_, a.real().vi_);
  EXPECT_FLOAT_EQ(5, a.real().val());
  EXPECT_FLOAT_EQ(-1, a.imag().val());
}
TEST(mathRevCore, stdComplexImag1) {
  cvar_t a(3, -1);
  EXPECT_FLOAT_EQ(-1, a.imag().val());
}
TEST(mathRevCore, stdComplexImag2) {
  // test double, int, var args
  cvar_t a(3, -1);
  a.imag(2.7);
  EXPECT_FLOAT_EQ(2.7, a.imag().val());
  EXPECT_FLOAT_EQ(3, a.real().val());

  a.imag(152);
  EXPECT_FLOAT_EQ(152, a.imag().val());
  EXPECT_FLOAT_EQ(3, a.real().val());

  var_t c(42);
  a.imag(c);
  EXPECT_FLOAT_EQ(42, a.imag().val());
  EXPECT_FLOAT_EQ(3, a.real().val());
}
TEST(mathRevCore, stdComplexOperatorEquals1and2and3) {
  auto f = [](auto& x, const auto& y) { x = y; };
  expect_compound_assign_operator(f);
}

TEST(mathRevCore, stdComplexOperatorAddEquals1and5) {
  auto f = [](auto& x1, const auto& x2) { return x1 += x2; };
  expect_compound_assign_operator(f);
}
TEST(mathRevCore, stdComplexOperatorSubtractEquals2and6) {
  auto f = [](auto& x1, const auto& x2) { return x1 -= x2; };
  expect_compound_assign_operator(f);
}
TEST(mathRevCore, stdComplexOperatorMultiplyEquals3and7) {
  auto f = [](auto& x1, const auto& x2) { return x1 *= x2; };
  expect_compound_assign_operator(f);
}
TEST(mathRevCore, stdComplexOperatorDivideEquals4and8) {
  auto f = [](auto& x1, const auto& x2) { return x1 /= x2; };
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
  EXPECT_FLOAT_EQ(bd, b.val());
}
TEST(mathRevCore, stdArgExternal1) {
  cdouble_t ad(1, 2);
  double bd = std::arg(ad);
  cvar_t a(1, 2);
  var_t b = std::arg(a);
  EXPECT_FLOAT_EQ(bd, b.val());
}
TEST(mathRevCore, stdNormExternal1) {
  cdouble_t ad(1, 2);
  double bd = std::norm(ad);
  cvar_t a(1, 2);
  var_t b = std::norm(a);
  EXPECT_FLOAT_EQ(bd, b.val());
}
TEST(mathRevCore, stdConj1) {
  cdouble_t ad(1, 2);
  cdouble_t bd = std::conj(ad);
  cvar_t a(1, 2);
  cvar_t b = std::conj(a);
  expect_complex(bd, b);
}
TEST(mathRevCore, stdProj1) {
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
TEST(mathRevCore, stdExp1) {
  expect_complex_common([](const auto& u) { return std::exp(u); });
}
TEST(mathRevCore, stdLog1) {
  expect_complex_common([](const auto& u) { return std::log(u); });
}
TEST(mathRevCore, stdLog101) {
  expect_complex_common([](const auto& u) { return std::log10(u); });
}
TEST(mathRevCore, stdPow1) {
  auto f = [](const auto& u, const auto& v) {
    // using std::pow;
    return pow(u, v);
  };
  expect_complex_common_binary(f);
  // verify (cvar_t, int) is unambiguous
  expect_complex(f(cdouble_t{1.2, 0.3}, 2), f(cvar_t{1.2, 0.3}, 2));
}
TEST(mathRevCore, stdSqrt1) {
  expect_complex_common([](const auto& u) { return std::sqrt(u); });
}
TEST(mathRevCore, stdSinh1) {
  expect_complex_common([](const auto& u) { return std::sinh(u); });
}
TEST(mathRevCore, stdCosh1) {
  expect_complex_common([](const auto& u) { return std::cosh(u); });
}
TEST(mathRevCore, stdTanh1) {
  expect_complex_common([](const auto& u) { return std::tanh(u); });
}
TEST(mathRevCore, stdAsinh1) {
  expect_complex_common([](const auto& u) { return std::asinh(u); });
}
TEST(mathRevCore, stdAcosh1) {
  expect_complex_common([](const auto& u) { return std::acosh(u); });
}
TEST(mathRevCore, stdAtanh1) {
  expect_complex_common([](const auto& u) { return std::atanh(u); });
}
TEST(mathRevCore, stdSin1) {
  expect_complex_common([](const auto& u) { return std::sin(u); });
}
TEST(mathRevCore, stdCos1) {
  expect_complex_common([](const auto& u) { return std::cos(u); });
}
TEST(mathRevCore, stdTan1) {
  expect_complex_common([](const auto& u) { return std::tan(u); });
}
TEST(mathRevCore, stdAsin1) {
  expect_complex_common([](const auto& u) { return std::asin(u); });
}
TEST(mathRevCore, stdAcos1) {
  expect_complex_common([](const auto& u) { return std::acos(u); });
}
TEST(mathRevCore, stdAtan1) {
  expect_complex_common([](const auto& u) { return std::acos(u); });
}

TEST(mathMix, iTimes) {
  auto f = [](const auto& x) {
    using stan::math::i_times;
    return to_array(i_times(from_array(x)));
  };
  auto zs = common_complex();
  for (const auto& z : zs) {
    stan::test::expect_ad(f, z);
  }
}
TEST(mathMix, negITimes) {
  auto f = [](const auto& x) {
    using stan::math::neg_i_times;
    return to_array(neg_i_times(from_array(x)));
  };
  auto zs = common_complex();
  for (const auto& z : zs) {
    stan::test::expect_ad(f, z);
  }
}

// type preserving/promoting to_complex
// adding 0.0 casts int to double but leaves others alone
template <typename U, typename V>
auto to_complex(const U& x, const V& y)
    -> std::complex<stan::return_type_t<U, V, double>> {
  return std::complex<stan::return_type_t<U, V, double>>(x, y);
}

template <typename U>
auto to_complex(const U& x) -> std::complex<stan::return_type_t<U, double>> {
  const auto y = x;
  return std::complex<stan::return_type_t<U, double>>(y);
}

TEST(mathMix, complexCtor1) {
  auto f = [](const auto& x) {
    auto a = to_complex(x);
    return to_array(a);
  };
  stan::test::expect_ad(f, 1.3);
  stan::test::expect_ad(f, 1);
}

TEST(mathMix, complexCtor2) {
  auto f = [](const auto& x, const auto& y) {
    auto a = to_complex(x, y);
    return to_array(a);
  };
  stan::test::expect_ad(f, 1.2, 1.3);
}
TEST(mathMix, operatorEqualScalar) {
  auto f = [](const auto& x) {
    std::complex<decltype(x + 0.0)> lhs;
    lhs = x;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2);
}
TEST(mathMix, operatorEqualComplex) {
  auto f = [](const auto& x, const auto& y) {
    std::complex<decltype(x + y + 0.0)> lhs(0.0, 0.0);
    auto rhs = to_complex(x, y);
    lhs = rhs;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, operatorPlusEqualScalar) {
  auto f = [](const auto& x) {
    std::complex<decltype(x + 0.0)> lhs;
    lhs += x;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2);
}
TEST(mathMix, operatorPlusEqualComplex) {
  auto f = [](const auto& x, const auto& y) {
    std::complex<decltype(x + y + 0.0)> lhs(1.5, -2.3);
    auto rhs = to_complex(x, y);
    lhs += rhs;
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, operatorMinusEqualScalar) {
  auto f = [](const auto& x) {
    std::complex<decltype(x + 0.0)> lhs;
    lhs -= x;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2);
}
TEST(mathMix, operatorMinusEqualComplex) {
  auto f = [](const auto& x, const auto& y) {
    std::complex<decltype(x + y + 0.0)> lhs(1.5, -2.3);
    auto rhs = to_complex(x, y);
    lhs -= rhs;
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, operatorTimesEqualScalar) {
  auto f = [](const auto& x) {
    std::complex<decltype(x + 0.0)> lhs;
    lhs *= x;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2);
}
TEST(mathMix, operatorTimesEqualComplex) {
  auto f = [](const auto& x, const auto& y) {
    std::complex<decltype(x + y + 0.0)> lhs(1.5, -2.3);
    auto rhs = to_complex(x, y);
    lhs *= rhs;
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, operatorDivideEqualScalar) {
  auto f = [](const auto& x) {
    std::complex<decltype(x + 0.0)> lhs;
    lhs /= x;  // assignment being tested
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2);
}
TEST(mathMix, operatorDivideEqualComplex) {
  auto f = [](const auto& x, const auto& y) {
    std::complex<decltype(x + y + 0.0)> lhs(1.5, -2.3);
    auto rhs = to_complex(x, y);
    lhs /= rhs;
    return to_array(lhs);
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, real) {
  auto f = [](const auto& x, const auto& y) {
    auto z = to_complex(x, y);
    return z.real();
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}
TEST(mathMix, imag) {
  auto f = [](const auto& x, const auto& y) {
    auto z = to_complex(x, y);
    return z.imag();
  };
  stan::test::expect_ad(f, 1.2, -2.17);
}

TEST(mathMix, operatorUnaryPlus) {
  auto f = [](const auto& x, const auto& y) {
    auto z = to_complex(x, y);
    return to_array(+z);
  };
  stan::test::expect_ad(f, 1.2, -3.1);
}
TEST(mathMix, operatorUnaryNegation) {
  auto f = [](const auto& x, const auto& y) {
    auto z = to_complex(x, y);
    return to_array(-z);
  };
  stan::test::expect_ad(f, 1.2, -3.1);
}

// TODO(carpenter): replace with param packs and move to test framework
std::vector<double> to_std_vec(double x1) { return {x1}; }
std::vector<double> to_std_vec(double x1, double x2) { return {x1, x2}; }
std::vector<double> to_std_vec(double x1, double x2, double x3) {
  return {x1, x2, x3};
}
std::vector<double> to_std_vec(double x1, double x2, double x3, double x4) {
  return {x1, x2, x3, x4};
}

// when all working refactor to this
// auto cwrap_zx = [](const auto& f) {
//   return [&](const auto& a, const auto & b) {
//     return to_array(f(from_array(a), b));
//   };
// };
// auto cwrap_xz = [](const auto& f) {
//   return [&](const auto& a, const auto & b) {
//     return to_array(f(a, from_array(b)));
//   };
// };
// auto cwrap_zz = [](const auto& f) {
//   return [&](const auto& a, const auto & b) {
//     return to_array(f(from_array(a), from_array(b)));
//   };
// };
// auto f = [](const auto& a, const auto& b) {
//   return a + b;
// };
// auto g = cwrap_zz(f);
// stan::test::expect_ad(g,
//                       to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

// FIX STARTING HERE!

TEST(mathMix, operatorAdd) {
  auto fzz = [](const auto& z1, const auto& z2) {
    auto c1 = from_array(z1);
    auto c2 = from_array(z2);
    auto y = c1 + c2;
    return to_array(y);
  };
  stan::test::expect_ad(fzz, to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

  auto fzx = [](const auto& z, const auto& x) {
    auto c = from_array(z);
    auto y = c + x;
    return to_array(y);
  };
  stan::test::expect_ad(fzx, to_std_vec(1.2, 2.3), -3.9);

  auto fxz = [](const auto& x, const auto& z) {
    auto c = from_array(z);
    auto y = x + c;
    return to_array(y);
  };
  stan::test::expect_ad(fxz, 1.2, to_std_vec(-3.9, -1.7));
}

TEST(mathMix, operatorSubtract) {
  auto fzz = [](const auto& z1, const auto& z2) {
    auto c1 = from_array(z1);
    auto c2 = from_array(z2);
    auto y = c1 - c2;
    return to_array(y);
  };
  stan::test::expect_ad(fzz, to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

  auto fzx = [](const auto& z, const auto& x) {
    auto c = from_array(z);
    auto y = c - x;
    return to_array(y);
  };
  stan::test::expect_ad(fzx, to_std_vec(1.2, 2.3), -3.9);

  auto fxz = [](const auto& x, const auto& z) {
    auto c = from_array(z);
    auto y = x - c;
    return to_array(y);
  };
  stan::test::expect_ad(fxz, 1.2, to_std_vec(-3.9, -1.7));
}

TEST(mathMix, operatorMultiply) {
  auto fzz = [](const auto& z1, const auto& z2) {
    auto c1 = from_array(z1);
    auto c2 = from_array(z2);
    auto y = c1 * c2;
    return to_array(y);
  };
  stan::test::expect_ad(fzz, to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

  auto fzx = [](const auto& z, const auto& x) {
    auto c = from_array(z);
    auto y = c * x;
    return to_array(y);
  };
  stan::test::expect_ad(fzx, to_std_vec(1.2, 2.3), -3.9);

  auto fxz = [](const auto& x, const auto& z) {
    auto c = from_array(z);
    auto y = x * c;
    return to_array(y);
  };
  stan::test::expect_ad(fxz, 1.2, to_std_vec(-3.9, -1.7));
}

TEST(mathMix, operatorDivide) {
  auto fzz = [](const auto& z1, const auto& z2) {
    auto c1 = from_array(z1);
    auto c2 = from_array(z2);
    auto y = c1 / c2;
    return to_array(y);
  };
  stan::test::expect_ad(fzz, to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

  auto fzx = [](const auto& z, const auto& x) {
    auto c = from_array(z);
    auto y = c / x;
    return to_array(y);
  };
  stan::test::expect_ad(fzx, to_std_vec(1.2, 2.3), -3.9);

  auto fxz = [](const auto& x, const auto& z) {
    auto c = from_array(z);
    auto y = x / c;
    return to_array(y);
  };
  stan::test::expect_ad(fxz, 1.2, to_std_vec(-3.9, -1.7));
}

TEST(mathMix, abs) {
  auto f = [](const auto& z) {
    using std::arg;
    return abs(from_array(z));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, arg) {
  auto f = [](const auto& z) {
    using std::arg;
    return arg(from_array(z));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, norm) {
  auto f = [](const auto& z) {
    using std::norm;
    return norm(from_array(z));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, conj) {
  auto f = [](const auto& z) {
    using std::norm;
    return to_array(conj(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, proj) {
  auto f = [](const auto& z) {
    using std::proj;
    return to_array(proj(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, exp) {
  auto f = [](const auto& z) {
    using std::exp;
    return to_array(exp(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, log) {
  auto f = [](const auto& z) {
    using std::log;
    return to_array(log(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, log10) {
  auto f = [](const auto& z) {
    using std::log10;
    return to_array(log10(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, pow) {
  auto fzz = [](const auto& z1, const auto& z2) {
    using std::pow;
    auto c1 = from_array(z1);
    auto c2 = from_array(z2);
    auto y = pow(c1, c2);
    return to_array(y);
  };
  stan::test::expect_ad(fzz, to_std_vec(1.2, 2.3), to_std_vec(-3.9, -1.7));

  auto fzx = [](const auto& z, const auto& x) {
    using std::pow;
    auto c = from_array(z);
    auto y = pow(c, x);
    return to_array(y);
  };
  stan::test::expect_ad(fzx, to_std_vec(1.2, 2.3), -3.9);

  auto fxz = [](const auto& x, const auto& z) {
    using std::pow;
    auto c = from_array(z);
    auto y = pow(x, c);
    return to_array(y);
  };
  stan::test::expect_ad(fxz, 1.2, to_std_vec(-3.9, -1.7));
}

TEST(mathMix, sqrt) {
  auto f = [](const auto& z) {
    using std::sqrt;
    return to_array(sqrt(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, sin) {
  auto f = [](const auto& z) {
    using std::sin;
    return to_array(sin(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, cos) {
  auto f = [](const auto& z) {
    using std::cos;
    return to_array(cos(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, tan) {
  auto f = [](const auto& z) {
    using std::tan;
    return to_array(tan(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, asin) {
  auto f = [](const auto& z) {
    using std::asin;
    return to_array(asin(from_array(z)));
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  stan::test::expect_ad(tols, f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, acos) {
  auto f = [](const auto& z) {
    using std::acos;
    return to_array(acos(from_array(z)));
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  stan::test::expect_ad(tols, f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, atan) {
  auto f = [](const auto& z) {
    using std::atan;
    return to_array(atan(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, sinh) {
  auto f = [](const auto& z) {
    using std::sinh;
    return to_array(sinh(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, cosh) {
  auto f = [](const auto& z) {
    using std::cosh;
    return to_array(cosh(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, tanh) {
  auto f = [](const auto& z) {
    using std::tanh;
    return to_array(tanh(from_array(z)));
  };
  stan::test::expect_ad(f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, asinh) {
  auto f = [](const auto& z) {
    using std::asinh;
    return to_array(asinh(from_array(z)));
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  stan::test::expect_ad(tols, f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, acosh) {
  auto f = [](const auto& z) {
    using std::acosh;
    return to_array(acosh(from_array(z)));
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  stan::test::expect_ad(tols, f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, atanh) {
  auto f = [](const auto& z) {
    using std::atanh;
    return to_array(atanh(from_array(z)));
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  stan::test::expect_ad(tols, f, to_std_vec(1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, -1.3));
  stan::test::expect_ad(tols, f, to_std_vec(-1.1, 2.9));
}

TEST(mathMix, polar) {
  auto f = [](const auto& r, const auto& theta) {
    using stan::math::complex_polar;
    return to_array(complex_polar(r, theta));
  };
  stan::test::expect_ad(f, 1.1, -1.3);
  stan::test::expect_ad(f, -1.1, -1.3);
  stan::test::expect_ad(f, -1.1, 2.9);
}

// tests of use, not of basic interface

std::vector<Eigen::MatrixXd> square_test_matrices() {
  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd a11(1, 1);
  a11 << -1.3;
  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 0.7;
  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 1.3, 0.7, 0.11, 0.13, -0.5, -1.7, 2.3;
  // a00 fails an assertion that we'd need to prevent
  //     if we implemented a function around this
  // a33 too hard numerically and adds 1m to tests
  return {a11, a22};
}
template <typename T>
void expectEigenSolver() {
  // test adapted from https://github.com/stan-dev/math/pull/789/files
  using matrix_v_t = Eigen::Matrix<T, -1, -1>;
  matrix_v_t a(3, 3);
  a << 1, 2, 3, 0.7, 0.11, 0.13, -5, -17, -23;
  Eigen::EigenSolver<matrix_v_t> s(a);
  auto ev = s.eigenvectors();
  auto I
      = (ev.inverse() * a * ev * s.eigenvalues().asDiagonal().inverse()).real();
  expect_identity_matrix(I);
}
template <typename T>
Eigen::EigenSolver<T> eigen_solver(const T& a) {
  Eigen::EigenSolver<T> s(a);
  return s;
}
TEST(mathMix, eigenSolver) {
  // value tests
  expectEigenSolver<var_t>();
  expectEigenSolver<fvar_d_t>();
  expectEigenSolver<fvar_fvar_d_t>();
  expectEigenSolver<fvar_v_t>();
  expectEigenSolver<fvar_fvar_v_t>();
  // derivative and value tests
  auto f1 = [](const auto& a) {
    return eigen_solver(a).eigenvectors().real().eval();
  };
  auto f2 = [](const auto& a) {
    return eigen_solver(a).eigenvectors().imag().eval();
  };
  auto g1 = [](const auto& a) {
    return eigen_solver(a).eigenvalues().real().eval();
  };
  auto g2 = [](const auto& a) {
    return eigen_solver(a).eigenvalues().imag().eval();
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 5e-3;
  tols.hessian_fvar_hessian_ = 5e-3;
  for (const auto& a : square_test_matrices()) {
    stan::test::expect_ad(tols, f1, a);
    stan::test::expect_ad(tols, f2, a);
    stan::test::expect_ad(tols, g1, a);
    stan::test::expect_ad(tols, g2, a);
  }
}

template <typename T>
void expectPseudoEigendecomposition() {
  // test adapted from https://github.com/stan-dev/math/pull/789/files
  using matrix_v_t = Eigen::Matrix<T, -1, -1>;
  matrix_v_t a(3, 3);
  a << 1, 2, 3, 0.7, 0.11, 0.13, -5, -17, -23;
  Eigen::EigenSolver<matrix_v_t> s(a);
  matrix_v_t D = s.pseudoEigenvalueMatrix();
  matrix_v_t V = s.pseudoEigenvectors();
  matrix_v_t I = V.inverse() * a * V * D.inverse();
  expect_identity_matrix(I);
}
TEST(mathMix, pseudoEigendecomposition) {
  expectPseudoEigendecomposition<var_t>();
  expectPseudoEigendecomposition<fvar_d_t>();
  expectPseudoEigendecomposition<fvar_fvar_d_t>();
  expectPseudoEigendecomposition<fvar_v_t>();
  expectPseudoEigendecomposition<fvar_fvar_v_t>();
  // the pseudo-eigendecomposition returns real-valued matrices
  auto f = [](const auto& a) {
    return eigen_solver(a).pseudoEigenvectors().eval();
  };
  auto g = [](const auto& a) {
    return eigen_solver(a).pseudoEigenvalueMatrix().eval();
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  for (const auto& a : square_test_matrices()) {
    stan::test::expect_ad(tols, f, a);
    stan::test::expect_ad(tols, g, a);
  }
}

template <typename T>
Eigen::ComplexSchur<T> complex_schur(const T& a) {
  Eigen::ComplexSchur<T> s(a);
  return s;
}
template <typename T>
void expectComplexSchur() {
  // test adapted from https://github.com/stan-dev/math/pull/789/files
  using matrix_v_t = Eigen::Matrix<T, -1, -1>;
  matrix_v_t a(3, 3);
  a << 1, 2, 3, 0.7, 0.11, 0.13, -5, -17, -23;
  Eigen::ComplexSchur<matrix_v_t> s(a);
  auto M = (s.matrixU().adjoint() * s.matrixU()).eval();
  matrix_v_t I = M.real() + M.imag();
  expect_identity_matrix(I);
}
TEST(mathMix, complexSchur) {
  expectComplexSchur<var_t>();
  expectComplexSchur<fvar_d_t>();
  expectComplexSchur<fvar_fvar_d_t>();
  expectComplexSchur<fvar_v_t>();
  expectComplexSchur<fvar_fvar_v_t>();
  auto f
      = [](const auto& a) { return complex_schur(a).matrixU().real().eval(); };
  auto g
      = [](const auto& a) { return complex_schur(a).matrixU().imag().eval(); };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  for (const auto& a : square_test_matrices()) {
    stan::test::expect_ad(tols, f, a);
    stan::test::expect_ad(tols, g, a);
  }
}

TEST(mathMix, powInt) {
  cvar_t v{3.0};
  cvar_t fv = pow(v, 2);
  var_t x{3.0};
  auto fx = pow(x, 2);
}

template <typename T>
void expect_complex_constructor() {
  std::complex<T> a{1};
  std::complex<T> b{2.0};
  std::complex<T> c{T(1.0)};
  std::complex<T> d{std::complex<double>(2.0)};
  std::complex<T> e{std::complex<T>(1.0)};
  SUCCEED();
}
TEST(mathMix, stdComplexConstructor) {
  using stan::math::fvar;
  using stan::math::var;
  expect_complex_constructor<double>();
  expect_complex_constructor<fvar<double>>();
  expect_complex_constructor<fvar<fvar<double>>>();
  expect_complex_constructor<var>();
  expect_complex_constructor<fvar<var>>();
  expect_complex_constructor<fvar<fvar<var>>>();
}

template <typename T>
void expect_complex_assignment() {
  std::complex<T> a = 1;
  std::complex<T> b = 2.0;
  std::complex<T> c = T(1.0);
  std::complex<T> d = std::complex<double>(2.0);
  std::complex<T> e = std::complex<T>(1.0);
  SUCCEED();
}
TEST(mathMix, stdComplexAssignment) {
  using stan::math::fvar;
  using stan::math::var;
  expect_complex_assignment<double>();
  expect_complex_assignment<fvar<double>>();
  expect_complex_assignment<fvar<fvar<double>>>();
  expect_complex_assignment<var>();
  expect_complex_assignment<fvar<var>>();
  expect_complex_assignment<fvar<fvar<var>>>();
}
// TEST(mathMix, traitsMeta) {
//   EXPECT_FALSE(stan::is_complex<int>::value);
//   EXPECT_FALSE(stan::is_complex<double>::value);
//   EXPECT_FALSE(stan::is_complex<stan::math::var>::value);
//   EXPECT_TRUE(stan::is_complex<std::complex<double>>::value);
//   EXPECT_TRUE(stan::is_complex<std::complex<stan::math::var>>::value);

//   EXPECT_TRUE(stan::is_arithmetic<int>::value);
//   EXPECT_TRUE(stan::is_arithmetic<double>::value);
//   EXPECT_TRUE(stan::is_arithmetic<stan::math::var>::value);
//   EXPECT_TRUE(stan::is_arithmetic<stan::math::fvar<double>>::value);
//   EXPECT_FALSE(stan::is_arithmetic<std::complex<double>>::value);
//   EXPECT_FALSE(stan::is_arithmetic<std::complex<stan::math::var>>::value);
//   EXPECT_FALSE(stan::is_arithmetic<std::string>::value);
// }

template <typename T>
void instantiate_multiply() {
  Eigen::Matrix<double, -1, -1> d(2, 2);
  d << 1, 2, 3, 4;
  Eigen::Matrix<T, -1, -1> v(2, 2);
  v << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<double>, -1, -1> cd(2, 2);
  cd << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<T>, -1, -1> cv(2, 2);
  cv << 1, 2, 3, 4;

  auto d_d = d * d;
  auto d_v = d * v;
  auto d_cd = d * cd;
  auto d_cv = d * cv;

  auto v_d = v * d;
  auto v_v = v * v;
  auto v_cd = v * cd;
  auto v_cv = v * cv;

  auto cd_d = cd * d;
  auto cd_v = cd * v;
  auto cd_cd = cd * cd;
  auto cd_cv = cd * cv;

  auto cv_d = cv * d;
  auto cv_v = cv * v;
  auto cv_cd = cv * cd;
  auto cv_cv = cv * cv;
}

TEST(mathMix, multiplicationPatterns) {
  using stan::math::fvar;
  using stan::math::var;
  instantiate_multiply<var>();
  instantiate_multiply<fvar<double>>();
  instantiate_multiply<fvar<fvar<double>>>();
  instantiate_multiply<fvar<var>>();
  instantiate_multiply<fvar<fvar<var>>>();
}

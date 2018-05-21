#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>
#include <limits>

TEST(AgradFwdAbs, Fvar) {
  using stan::math::fvar;
  using std::abs;
  using std::isnan;

  fvar<double> x(2);
  fvar<double> y(-3);
  x.d_ = 1.0;
  y.d_ = 2.0;

  fvar<double> a = abs(x);
  EXPECT_FLOAT_EQ(abs(2), a.val_);
  EXPECT_FLOAT_EQ(1.0, a.d_);

  fvar<double> b = abs(-x);
  EXPECT_FLOAT_EQ(abs(-2), b.val_);
  EXPECT_FLOAT_EQ(1.0, b.d_);

  fvar<double> c = abs(y);
  EXPECT_FLOAT_EQ(abs(-3), c.val_);
  EXPECT_FLOAT_EQ(-2.0, c.d_);

  fvar<double> d = abs(2 * x);
  EXPECT_FLOAT_EQ(abs(2 * 2), d.val_);
  EXPECT_FLOAT_EQ(2 * 1.0, d.d_);

  fvar<double> e = abs(y + 4);
  EXPECT_FLOAT_EQ(abs(-3 + 4), e.val_);
  EXPECT_FLOAT_EQ(2.0, e.d_);

  fvar<double> f = abs(x - 2);
  EXPECT_FLOAT_EQ(abs(2 - 2), f.val_);
  EXPECT_FLOAT_EQ(0, f.d_);

  fvar<double> z = std::numeric_limits<double>::quiet_NaN();
  fvar<double> g = abs(z);
  EXPECT_TRUE(boost::math::isnan(g.val_));
  EXPECT_TRUE(boost::math::isnan(g.d_));

  fvar<double> w = 0;
  fvar<double> h = abs(w);
  EXPECT_FLOAT_EQ(0.0, h.val_);
  EXPECT_FLOAT_EQ(0.0, h.d_);
}

TEST(AgradFwdAbs, complexNotNullIssue123) {
  using stan::math::fvar;
  fvar<double> x(3, 1);
  std::complex<fvar<double>> z = x;
  EXPECT_TRUE(z.real().val());
  EXPECT_FALSE(z.imag().val());
  auto y(z.real() + 1);
  // gradient propagates partly through imaginary <- y <- real <- x
  z.imag(y);
  EXPECT_TRUE(z.imag().val());
  auto zabs(abs(z));
  EXPECT_EQ(zabs.val(), 5);
  // z^2 = x^2 + (x+1)^2 = 2x^2+2x+1 ==> 2z dz/dx = 4x + 2
  // dz/dx = (4x + 2)/(2 z) = (4 * 3 + 2)/(2 * 5) = 1.4
  EXPECT_EQ(zabs.tangent(), 1.4);

  // the following are complex var left and right op parameter coercions
  // (i.e. checks for correct instantiations)

  // 8/((1+3+4i)*2) = 1/(1+1i) * (1-i)/(1-i) = (1-i)/2
  auto q(8 / ((1 + z) * 2));
  q += std::complex<double>(.5, 1.5);             // 0.5-0.5i + .5+1.5i = 1+1i
  std::complex<fvar<double>> r(2 * (z + 1) / 8);  // 2*(3+4i+1)/8  = 1+1i
  EXPECT_TRUE(q == std::complex<double>(1, 1));
  EXPECT_TRUE(std::complex<double>(1, 1) == q);
  EXPECT_TRUE(r == std::complex<double>(1, 1));
  EXPECT_TRUE(std::complex<double>(1, 1) == r);
  EXPECT_TRUE(q == r);
  EXPECT_TRUE(r == q);

  EXPECT_NEAR((abs(q) * abs(r)).val(), 2, 1E-8);
  EXPECT_NEAR((abs(r) * abs(q)).val(), 2, 1E-8);
}

TEST(AgradFwdAbs, FvarFvarDouble) {
  using stan::math::fvar;
  using std::abs;

  fvar<fvar<double>> x;
  x.val_.val_ = 4.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double>> a = abs(x);

  EXPECT_FLOAT_EQ(4.0, a.val_.val_);
  EXPECT_FLOAT_EQ(1.0, a.val_.d_);
  EXPECT_FLOAT_EQ(0.0, a.d_.val_);
  EXPECT_FLOAT_EQ(0.0, a.d_.d_);
}

struct abs_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return abs(arg1);
  }
};

TEST(AgradFwdAbs, abs_NaN) {
  abs_fun abs_;
  test_nan_fwd(abs_, false);
}

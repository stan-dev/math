#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>
#include <limits>

TEST(MathFunctions, inv_Phi) {
  using stan::math::Phi;
  using stan::math::fvar;
  using stan::math::inv_Phi;
  EXPECT_FLOAT_EQ(0.0, inv_Phi(0.5));
  fvar<double> p = 0.123456789;
  EXPECT_FLOAT_EQ(p.val_, Phi(inv_Phi(p)).val_);
  p = 8e-311;
  EXPECT_FLOAT_EQ(p.val_, Phi(inv_Phi(p)).val_);
  p = 0.99;
  EXPECT_FLOAT_EQ(p.val_, Phi(inv_Phi(p)).val_);

  // breakpoints
  p = 0.02425;
  EXPECT_FLOAT_EQ(p.val_, Phi(inv_Phi(p)).val_);
  p = 0.97575;
  EXPECT_FLOAT_EQ(p.val_, Phi(inv_Phi(p)).val_);
}
TEST(MathFunctions, inv_Phi_inf) {
  using stan::math::fvar;
  using stan::math::inv_Phi;
  fvar<double> p = 7e-311;
  const fvar<double> inf = std::numeric_limits<fvar<double> >::infinity();
  EXPECT_EQ(inv_Phi(p), -inf);
  p = 1.0;
  EXPECT_EQ(inv_Phi(p), inf);
}
TEST(MathFunctions, inv_Phi_nan) {
  using stan::math::fvar;
  using stan::math::inv_Phi;
  fvar<double> nan = std::numeric_limits<fvar<double> >::quiet_NaN();
  EXPECT_THROW(inv_Phi(nan), std::domain_error);
  EXPECT_THROW(inv_Phi(-2.0), std::domain_error);
  EXPECT_THROW(inv_Phi(2.0), std::domain_error);
}

TEST(AgradFwdinv_Phi, Fvar) {
  using stan::math::Phi;
  using stan::math::fvar;
  using stan::math::inv_Phi;

  fvar<double> x = 0.1;
  x.d_ = 1.0;

  fvar<double> y = Phi(inv_Phi(x));

  EXPECT_FLOAT_EQ(0.1, y.val_);
  EXPECT_FLOAT_EQ(1.0, y.d_);
}
TEST(AgradFwdinv_Phi, FvarFvarDouble) {
  using stan::math::Phi;
  using stan::math::fvar;
  using stan::math::inv_Phi;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.1;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = Phi(inv_Phi(x));

  EXPECT_FLOAT_EQ(0.1, a.val_.val_);
  EXPECT_FLOAT_EQ(1.0, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct inv_Phi_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return inv_Phi(arg1);
  }
};

TEST(AgradFwdinv_Phi, inv_Phi_NaN) {
  inv_Phi_fun foo;
  test_nan_fwd(foo, true);
}

#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>
#include <stdexcept>

TEST(AgradFwdLog1m,Fvar) {
  using stan::math::fvar;
  using stan::math::log1m;
  using std::isnan;

  fvar<double> x(0.5,1.0);

  fvar<double> a = log1m(x);
  EXPECT_FLOAT_EQ(log1m(0.5), a.val_);
  EXPECT_FLOAT_EQ(-1 / (1 - 0.5), a.d_);
}
TEST(AgradFwdLog1m, FvarExcepts) {
  using stan::math::fvar;
  using stan::math::log1m;
  EXPECT_THROW(stan::math::log1m(fvar<double>(2)), std::domain_error);
}

TEST(MathFunctions, log1m_inf_return) {
  using stan::math::fvar;
  using stan::math::log1m;
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            log1m(fvar<double>(1.0)));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            log1m(fvar<double>(1)));
}

TEST(AgradFwdLog1m,FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::log1m;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = log1m(x);

  EXPECT_FLOAT_EQ(log1m(0.5), a.val_.val_);
  EXPECT_FLOAT_EQ(-1 / (0.5), a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  a = log1m(y);
  EXPECT_FLOAT_EQ(log1m(0.5), a.val_.val_);
  EXPECT_FLOAT_EQ(0, a.val_.d_);
  EXPECT_FLOAT_EQ(-1 / (0.5), a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct log1m_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log1m(arg1);
  }
};

TEST(AgradFwdLog1m,log1m_NaN) {
  log1m_fun log1m_;
  test_nan_fwd(log1m_,false);
}

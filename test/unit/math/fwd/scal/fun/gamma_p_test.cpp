#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdGammaP, gamma_p) {
  using boost::math::gamma_p;
  using stan::math::fvar;
  using stan::math::gamma_p;

  fvar<double> x(0.5001);
  x.d_ = 1.0;
  fvar<double> y(1.0001);
  y.d_ = 1.0;

  fvar<double> a = gamma_p(x, y);
  EXPECT_FLOAT_EQ(gamma_p(0.5001, 1.0001), a.val_);
  EXPECT_FLOAT_EQ(
      boost::math::gamma_p_derivative(0.5001, 1.0001) - 0.3898178624664172,
      a.d_);

  double z = 1.0001;
  double w = 0.5001;

  a = gamma_p(x, z);
  EXPECT_FLOAT_EQ(gamma_p(0.5001, 1.0001), a.val_);
  EXPECT_FLOAT_EQ(-0.3898178624664172, a.d_);

  a = gamma_p(w, y);
  EXPECT_FLOAT_EQ(gamma_p(0.5001, 1.0001), a.val_);
  EXPECT_FLOAT_EQ(boost::math::gamma_p_derivative(0.5001, 1.0001), a.d_);

  EXPECT_THROW(gamma_p(-x, y), std::domain_error);
  EXPECT_THROW(gamma_p(x, -y), std::domain_error);
}

TEST(AgradFwdGammaP, FvarFvarDouble) {
  using boost::math::gamma_p;
  using stan::math::fvar;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.5001;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 1.0001;
  y.d_.val_ = 1.0;

  fvar<fvar<double> > a = gamma_p(x, y);

  EXPECT_FLOAT_EQ(gamma_p(0.5001, 1.0001), a.val_.val_);
  EXPECT_FLOAT_EQ(-0.3898178624664172, a.val_.d_);
  EXPECT_FLOAT_EQ(boost::math::gamma_p_derivative(0.5001, 1.0001), a.d_.val_);

  EXPECT_FLOAT_EQ(0.40747109, a.d_.d_);
}

struct gamma_p_fun {
  template <typename T0, typename T1>
  inline typename boost::math::tools::promote_args<T0, T1>::type operator()(
      const T0 arg1, const T1 arg2) const {
    return gamma_p(arg1, arg2);
  }
};

TEST(AgradFwdGammaP, nan) {
  gamma_p_fun gamma_p_;
  test_nan_fwd(gamma_p_, 3.0, 5.0, false);
}

#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/digamma.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdBinomialCoefficientLog,Fvar) {
  using stan::math::fvar;
  using stan::math::binomial_coefficient_log;
  using boost::math::digamma;

  fvar<double> x(2004.0,1.0);
  fvar<double> y(1002.0,2.0);

  fvar<double> a = stan::math::binomial_coefficient_log(x, y);
  EXPECT_FLOAT_EQ(binomial_coefficient_log(2004.0, 1002.0), a.val_);
  EXPECT_FLOAT_EQ(0.69289774, a.d_);
}


TEST(AgradFwdBinomialCoefficientLog,FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::binomial_coefficient_log;
  using stan::math::binomial_coefficient_log;

  fvar<fvar<double> > x;
  x.val_.val_ = 2004.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 1002.0;
  y.d_.val_ = 1.0;

  fvar<fvar<double> > a = binomial_coefficient_log(x,y);

  EXPECT_FLOAT_EQ(binomial_coefficient_log(2004.0,1002.0), a.val_.val_);
  EXPECT_FLOAT_EQ(0.69289774, a.val_.d_);
  EXPECT_NEAR(0, a.d_.val_,1e-8);
  EXPECT_FLOAT_EQ(0.0009975062, a.d_.d_);
}


struct binomial_coefficient_log_fun {
  template <typename T0, typename T1>
  inline 
  typename boost::math::tools::promote_args<T0,T1>::type
  operator()(const T0 arg1,
             const T1 arg2) const {
    return binomial_coefficient_log(arg1,arg2);
  }
};

TEST(AgradFwdBinomialCoefficientLog, nan) {
  binomial_coefficient_log_fun binomial_coefficient_log_;
  test_nan_fwd(binomial_coefficient_log_,3.0,5.0,false);
}

#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/digamma.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdLogFallingFactorial,Fvar) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using boost::math::digamma;

  fvar<double> a(4.0, 1.0);
  fvar<double> x = log_falling_factorial(a, 2);
  EXPECT_FLOAT_EQ(std::log(12.0), x.val_);
  EXPECT_FLOAT_EQ((boost::math::digamma(5) - boost::math::digamma(3)),
                  x.d_);
  
  //finite diff
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(4.0 + eps, 2.0)
                  - stan::math::log_falling_factorial(4.0 - eps, 2.0))
                  / (2 * eps), x.d_);

  fvar<double> c(-3.0, 2.0);

  EXPECT_THROW(log_falling_factorial(c, 2), std::domain_error);
  EXPECT_THROW(log_falling_factorial(c, c), std::domain_error);

  x = log_falling_factorial(a,a);
  EXPECT_FLOAT_EQ(std::log(24.0), x.val_);
  EXPECT_FLOAT_EQ(boost::math::digamma(5), x.d_);

  x = log_falling_factorial(5, a);
  EXPECT_FLOAT_EQ(std::log(120.0), x.val_);
  EXPECT_FLOAT_EQ(digamma(2.0),x.d_);
  
  //finite diff
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(5.0, 4.0 + eps)
                  - stan::math::log_falling_factorial(5.0, 4.0 - eps))
                  / (2 * eps), x.d_);
}

TEST(AgradFwdLogFallingFactorial,FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;

  fvar<fvar<double> > x;
  x.val_.val_ = 4.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 3.0;
  y.d_.val_ = 1.0;

  fvar<fvar<double> > a = log_falling_factorial(x,y);

  EXPECT_FLOAT_EQ(3.1780539, a.val_.val_);
  EXPECT_FLOAT_EQ(1.0833334, a.val_.d_);
  EXPECT_FLOAT_EQ(0.42278433, a.d_.val_);
  EXPECT_FLOAT_EQ(0.64493406, a.d_.d_);
}

struct log_falling_factorial_fun {
  template <typename T0, typename T1>
  inline 
  typename boost::math::tools::promote_args<T0,T1>::type
  operator()(const T0 arg1,
             const T1 arg2) const {
    return log_falling_factorial(arg1,arg2);
  }
};

TEST(AgradFwdLogFallingFactorial, nan) {
  log_falling_factorial_fun log_falling_factorial_;
  test_nan_fwd(log_falling_factorial_,3.0,5.0,false);
}

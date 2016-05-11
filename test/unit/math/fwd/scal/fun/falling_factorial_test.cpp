#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdFallingFactorial,Fvar) {
  using stan::math::fvar;
  using stan::math::falling_factorial;
  using boost::math::digamma;

  fvar<double> a(4.0, 1.0);
  fvar<double> x = falling_factorial(a, 2);
  EXPECT_FLOAT_EQ(12.0, x.val_);
  EXPECT_FLOAT_EQ((boost::math::digamma(5) - boost::math::digamma(3))
                  * 12.0, x.d_);
  
  //finite diff
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::falling_factorial(4.0 + eps, 2.0)
                  - stan::math::falling_factorial(4.0 - eps, 2.0))
                  / (2 * eps), x.d_);

  fvar<double> c(-3.0, 2.0);

  EXPECT_THROW(falling_factorial(c, 2), std::domain_error);
  EXPECT_THROW(falling_factorial(c, c), std::domain_error);

  x = falling_factorial(a,a);
  EXPECT_FLOAT_EQ(24.0, x.val_);
  EXPECT_FLOAT_EQ(24.0 * boost::math::digamma(5), x.d_);

  x = falling_factorial(5, a);
  EXPECT_FLOAT_EQ(120.0, x.val_);
  EXPECT_FLOAT_EQ(120.0 * digamma(2.0),x.d_);
  
  //finite diff
  EXPECT_FLOAT_EQ((stan::math::falling_factorial(5.0, 4.0 + eps)
                  - stan::math::falling_factorial(5.0, 4.0 - eps))
                  / (2 * eps), x.d_);
}

TEST(AgradFwdFallingFactorial,FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::falling_factorial;

  fvar<fvar<double> > x;
  x.val_.val_ = 4.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 4.0;
  y.d_.val_ = 1.0;

  fvar<fvar<double> > a = falling_factorial(x,y);

  EXPECT_FLOAT_EQ(falling_factorial(4,4.0), a.val_.val_);
  EXPECT_FLOAT_EQ(50, a.val_.d_);
  EXPECT_FLOAT_EQ(-13.853176, a.d_.val_);
  EXPECT_FLOAT_EQ(10.617635, a.d_.d_);
}
struct falling_factorial_fun {
  template <typename T0, typename T1>
  inline 
  typename boost::math::tools::promote_args<T0,T1>::type
  operator()(const T0 arg1,
             const T1 arg2) const {
    return falling_factorial(arg1,arg2);
  }
};

TEST(AgradFwdFallingFactorial, nan) {
  falling_factorial_fun falling_factorial_;
  test_nan_fwd(falling_factorial_,3.0,5.0,false);
}

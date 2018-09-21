#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdRisingFactorial, Fvar) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;

  fvar<double> a(4.0, 1.0);
  fvar<double> x = rising_factorial(a, 1);
  EXPECT_FLOAT_EQ(4.0, x.val_);
  EXPECT_FLOAT_EQ(1.0, x.d_);

  // finite diff
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(4.0 + eps, 1)
                   - stan::math::rising_factorial(4.0 - eps, 1))
                      / (2 * eps),
                  x.d_);

  fvar<double> c(-3.0, 2.0);

  EXPECT_THROW(rising_factorial(c, -2), std::domain_error);
}

TEST(AgradFwdRisingFactorial, FvarFvarDouble) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;

  fvar<fvar<double> > x;
  x.val_.val_ = 4.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = rising_factorial(x, 4);

  EXPECT_FLOAT_EQ((840.0), a.val_.val_);
  EXPECT_FLOAT_EQ(840. * (digamma(8) - digamma(4)), a.val_.d_);
}

struct rising_factorial_fun {
  template <typename T>
  inline typename boost::math::tools::promote_args<T>::type operator()(
      const T arg1, int arg2) const {
    return rising_factorial(arg1, arg2);
  }
};

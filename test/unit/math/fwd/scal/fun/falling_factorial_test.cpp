#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdFallingFactorial, Fvar) {
  using stan::math::digamma;
  using stan::math::falling_factorial;
  using stan::math::fvar;

  fvar<double> a(4.0, 1.0);
  fvar<double> x = falling_factorial(a, 2);
  EXPECT_FLOAT_EQ(12.0, x.val_);
  EXPECT_FLOAT_EQ((stan::math::digamma(5) - stan::math::digamma(3)) * 12.0,
                  x.d_);

  // finite diff
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::falling_factorial(4.0 + eps, 2)
                   - stan::math::falling_factorial(4.0 - eps, 2))
                      / (2 * eps),
                  x.d_);

  fvar<double> c(-3.0, 2.0);

  EXPECT_THROW(falling_factorial(c, -2), std::domain_error);
}

TEST(AgradFwdFallingFactorial, FvarFvarDouble) {
  using stan::math::falling_factorial;
  using stan::math::fvar;

  fvar<fvar<double> > x;
  x.val_.val_ = 4.0;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = falling_factorial(x, 4);

  EXPECT_FLOAT_EQ(falling_factorial(4, 4), a.val_.val_);
  EXPECT_FLOAT_EQ(50, a.val_.d_);
}
struct falling_factorial_fun {
  template <typename T>
  inline typename boost::math::tools::promote_args<T>::type operator()(
      const T arg1, int arg2) const {
    return falling_factorial(arg1, arg2);
  }
};

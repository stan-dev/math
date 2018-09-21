#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/digamma.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

TEST(AgradFwdFallingFactorial, FvarVar_1stDeriv) {
  using stan::math::digamma;
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a(5.0, 1.0);
  fvar<var> c = falling_factorial(a, 3);

  EXPECT_FLOAT_EQ(falling_factorial(5, 3), c.val_.val());
  EXPECT_FLOAT_EQ(
      falling_factorial(5, 3) * (digamma(5 + 1) - digamma(5 + 1 - 3)),
      c.d_.val());

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.val_.grad(y, g);
  EXPECT_FLOAT_EQ(47, g[0]);
  EXPECT_FLOAT_EQ(0, g[1]);
}

TEST(AgradFwdFallingFactorial, FvarVar_2ndDeriv_x) {
  using stan::math::digamma;
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::trigamma;
  using stan::math::var;
  using std::pow;

  fvar<var> a(5.0, 1.0);
  fvar<var> c = falling_factorial(a, 3);

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(falling_factorial(5, 3)
                  * (pow((digamma(5 + 1) - digamma(5 + 1 - 3)), 2)
                     - trigamma(5 + 1 - 3) + trigamma(5 + 1)),
              g[0], 0.1);
}

TEST(AgradFwdFallingFactorial, FvarVar_2ndDeriv_y) {
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::var;

  /**
   * Second derivative w.r.t. n should return 0,
   * since n is an integer.
   */

  fvar<var> a(5.0, 0.0);
  fvar<var> c = falling_factorial(a, 3);

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(0, g[1], 0.1);
}

TEST(AgradFwdFallingFactorial, FvarFvarVar_1stDeriv) {
  using stan::math::digamma;
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = falling_factorial(x, 3);

  EXPECT_FLOAT_EQ(falling_factorial(5, 3), a.val_.val_.val());
  EXPECT_FLOAT_EQ(47, a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  ASSERT_NEAR(0, a.d_.d_.val(), .01);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.val_.val_.grad(p, g);
  EXPECT_FLOAT_EQ(
      falling_factorial(5, 3) * (digamma(5 + 1) - digamma(5 + 1 - 3)), g[0]);
  EXPECT_FLOAT_EQ(0, g[1]);
}

TEST(AgradFwdFallingFactorial, FvarFvarVar_2ndDeriv_x) {
  using stan::math::digamma;
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::trigamma;
  using stan::math::var;
  using std::pow;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = falling_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.val_.d_.grad(p, g);
  ASSERT_NEAR(falling_factorial(5, 3)
                  * (pow((digamma(5 + 1) - digamma(5 + 1 - 3)), 2)
                     - trigamma(5 + 1 - 3) + trigamma(5 + 1)),
              g[0], 0.01);
  ASSERT_NEAR(0, g[1], 0.01);
}

TEST(AgradFwdFallingFactorial, FvarFvarVar_2ndDeriv_y) {
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = falling_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.d_.val_.grad(p, g);
  ASSERT_NEAR(0, g[0], 0.01);
  ASSERT_NEAR(0, g[1], 0.01);
}

TEST(AgradFwdFallingFactorial, FvarFvarVar_3rdDeriv) {
  using stan::math::falling_factorial;
  using stan::math::fvar;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = falling_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.d_.d_.grad(p, g);
  ASSERT_NEAR(0, g[0], 0.03);
  ASSERT_NEAR(0, g[1], 0.03);
}

struct falling_factorial_fun {
  template <typename T>
  inline typename boost::math::tools::promote_args<T>::type operator()(
      const T arg1, int arg2) const {
    return falling_factorial(arg1, arg2);
  }
};

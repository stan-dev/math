#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

TEST(AgradFwdRisingFactorial, FvarVar_1stDeriv) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::var;

  fvar<var> a(5.0, 1.0);
  fvar<var> c = rising_factorial(a, 3);

  EXPECT_FLOAT_EQ(rising_factorial(5, 3), c.val_.val());
  EXPECT_FLOAT_EQ(rising_factorial(5, 3) * (digamma(5 + 3) - digamma(5)),
                  c.d_.val());

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.val_.grad(y, g);
  EXPECT_FLOAT_EQ(rising_factorial(5, 3) * (digamma(5 + 3) - digamma(5)), g[0]);
  EXPECT_FLOAT_EQ(0, g[1]);
}

TEST(AgradFwdRisingFactorial, FvarVar_2ndDeriv_x) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::trigamma;
  using stan::math::var;
  using std::pow;

  fvar<var> a(5.0, 1.0);
  fvar<var> c = rising_factorial(a, 3);

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(rising_factorial(5, 3)
                  * (pow((digamma(5 + 3) - digamma(5)), 2) + trigamma(5 + 3)
                     - trigamma(5)),
              g[0], 0.1);
}

TEST(AgradFwdRisingFactorial, FvarVar_2ndDeriv_y) {
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::var;

  /**
   * Second derivative w.r.t. n should return 0,
   * since n is an integer.
   */

  fvar<var> a(5.0, 0.0);
  fvar<var> c = rising_factorial(a, 3);

  AVEC y = createAVEC(a.val_, 3);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(0, g[1], 0.1);
}

TEST(AgradFwdRisingFactorial, FvarFvarVar_1stDeriv) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = rising_factorial(x, 3);

  EXPECT_FLOAT_EQ(rising_factorial(5, 3), a.val_.val_.val());
  EXPECT_FLOAT_EQ(107.0, a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  ASSERT_NEAR(0, a.d_.d_.val(), .01);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.val_.val_.grad(p, g);
  EXPECT_FLOAT_EQ(rising_factorial(5, 3) * (digamma(5 + 3) - digamma(5)), g[0]);
  EXPECT_FLOAT_EQ(0, g[1]);
}

TEST(AgradFwdRisingFactorial, FvarFvarVar_2ndDeriv_x) {
  using stan::math::digamma;
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::trigamma;
  using stan::math::var;
  using std::pow;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = rising_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.val_.d_.grad(p, g);
  ASSERT_NEAR(rising_factorial(5, 3)
                  * (pow((digamma(5 + 3) - digamma(5)), 2) + trigamma(5 + 3)
                     - trigamma(5)),
              g[0], 0.01);
  ASSERT_NEAR(0, g[1], 0.01);
}

TEST(AgradFwdRisingFactorial, FvarFvarVar_2ndDeriv_y) {
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = rising_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.d_.val_.grad(p, g);
  ASSERT_NEAR(0, g[0], 0.01);
  ASSERT_NEAR(0, g[1], 0.01);
}

TEST(AgradFwdRisingFactorial, FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::rising_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = rising_factorial(x, 3);

  AVEC p = createAVEC(x.val_.val_, 3);
  VEC g;
  a.d_.d_.grad(p, g);
  ASSERT_NEAR(0, g[0], 0.03);
  ASSERT_NEAR(0, g[1], 0.03);
}

struct rising_factorial_fun {
  template <typename T>
  inline typename boost::math::tools::promote_args<T>::type operator()(
      const T arg1, int arg2) const {
    return rising_factorial(arg1, arg2);
  }
};

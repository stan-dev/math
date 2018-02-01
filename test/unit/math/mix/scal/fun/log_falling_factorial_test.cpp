#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/digamma.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

double eps = 1e-6;
double first_deriv_a = (stan::math::log_falling_factorial(5.0 + eps, 3.0)
                        - stan::math::log_falling_factorial(5.0 - eps, 3.0))
                       / (2 * eps);
double first_deriv_b = (stan::math::log_falling_factorial(5.0, 3.0 + eps)
                        - stan::math::log_falling_factorial(5.0, 3.0 - eps))
                       / (2 * eps);

double eps2 = 1e-4;
double second_deriv_aa
    = (stan::math::log_falling_factorial(5.0 + 2 * eps2, 3.0)
       - 2 * stan::math::log_falling_factorial(5.0 + eps2, 3.0)
       + stan::math::log_falling_factorial(5.0, 3.0))
      / std::pow(eps2, 2);
double second_deriv_bb
    = (stan::math::log_falling_factorial(5.0, 3.0 + 2 * eps2)
       - 2 * stan::math::log_falling_factorial(5.0, 3.0 + eps2)
       + stan::math::log_falling_factorial(5.0, 3.0))
      / std::pow(eps2, 2);
double second_deriv_ab
    = (stan::math::log_falling_factorial(5.0 + eps2, 3.0 + eps2)
       - stan::math::log_falling_factorial(5.0 - eps2, 3.0 + eps2)
       - stan::math::log_falling_factorial(5.0 + eps2, 3.0 - eps2)
       + stan::math::log_falling_factorial(5.0 - eps2, 3.0 - eps2))
      / 4 / std::pow(eps2, 2);

double third_deriv_aab
    = (stan::math::log_falling_factorial(5.0 + 2 * eps2, 3.0 + eps2)
       - 2 * stan::math::log_falling_factorial(5.0 + eps2, 3.0 + eps2)
       + stan::math::log_falling_factorial(5.0, 3.0 + eps2)
       - stan::math::log_falling_factorial(5.0 + 2 * eps2, 3.0 - eps2)
       + 2 * stan::math::log_falling_factorial(5.0 + eps2, 3.0 - eps2)
       - stan::math::log_falling_factorial(5.0, 3.0 - eps2))
      / 2 / std::pow(eps2, 3);

double third_deriv_abb
    = (stan::math::log_falling_factorial(5.0 + eps2, 3.0 + 2 * eps2)
       - 2 * stan::math::log_falling_factorial(5.0 + eps2, 3.0 + eps2)
       + stan::math::log_falling_factorial(5.0 + eps2, 3.0)
       - stan::math::log_falling_factorial(5.0 - eps2, 3.0 + 2 * eps2)
       + 2 * stan::math::log_falling_factorial(5.0 - eps2, 3.0 + eps2)
       - stan::math::log_falling_factorial(5.0 - eps2, 3.0))
      / 2 / std::pow(eps2, 3);

TEST(AgradFwdLogFallingFactorial, FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<var> a(5.0, 1.0);
  fvar<var> b(3.0, 1.0);
  fvar<var> c = log_falling_factorial(a, b);

  EXPECT_FLOAT_EQ(std::log(60), c.val_.val());
  EXPECT_FLOAT_EQ(first_deriv_a + first_deriv_b, c.d_.val());

  AVEC y = createAVEC(a.val_, b.val_);
  VEC g;
  c.val_.grad(y, g);
  EXPECT_FLOAT_EQ(first_deriv_a, g[0]);
  EXPECT_FLOAT_EQ(first_deriv_b, g[1]);
}
TEST(AgradFwdLogFallingFactorial, FvarVar_2ndDeriv_x) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<var> a(5.0, 1.0);
  fvar<var> b(3.0, 0.0);
  fvar<var> c = log_falling_factorial(a, b);

  AVEC y = createAVEC(a.val_, b.val_);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(second_deriv_aa, g[0], 0.1);
}
TEST(AgradFwdLogFallingFactorial, FvarVar_2ndDeriv_y) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<var> a(5.0, 0.0);
  fvar<var> b(3.0, 1.0);
  fvar<var> c = log_falling_factorial(a, b);

  AVEC y = createAVEC(a.val_, b.val_);
  VEC g;
  c.d_.grad(y, g);
  ASSERT_NEAR(second_deriv_bb, g[1], 0.1);
}
TEST(AgradFwdLogFallingFactorial, FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > y;
  y.val_.val_ = 3.0;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = log_falling_factorial(x, y);

  EXPECT_FLOAT_EQ(log_falling_factorial(5, 3.0), a.val_.val_.val());
  EXPECT_FLOAT_EQ(first_deriv_a, a.val_.d_.val());
  EXPECT_FLOAT_EQ(first_deriv_b, a.d_.val_.val());
  ASSERT_NEAR(second_deriv_ab, a.d_.d_.val(), .01);

  AVEC p = createAVEC(x.val_.val_, y.val_.val_);
  VEC g;
  a.val_.val_.grad(p, g);
  EXPECT_FLOAT_EQ(first_deriv_a, g[0]);
  EXPECT_FLOAT_EQ(first_deriv_b, g[1]);
}

TEST(AgradFwdLogFallingFactorial, FvarFvarVar_2ndDeriv_x) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > y;
  y.val_.val_ = 3.0;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = log_falling_factorial(x, y);

  AVEC p = createAVEC(x.val_.val_, y.val_.val_);
  VEC g;
  a.val_.d_.grad(p, g);
  ASSERT_NEAR(second_deriv_aa, g[0], 0.01);
  ASSERT_NEAR(second_deriv_ab, g[1], 0.01);
}
TEST(AgradFwdLogFallingFactorial, FvarFvarVar_2ndDeriv_y) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > y;
  y.val_.val_ = 3.0;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = log_falling_factorial(x, y);

  AVEC p = createAVEC(x.val_.val_, y.val_.val_);
  VEC g;
  a.d_.val_.grad(p, g);
  ASSERT_NEAR(second_deriv_ab, g[0], 0.01);
  ASSERT_NEAR(second_deriv_bb, g[1], 0.01);
}
TEST(AgradFwdLogFallingFactorial, FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::log_falling_factorial;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 5.0;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > y;
  y.val_.val_ = 3.0;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = log_falling_factorial(x, y);

  AVEC p = createAVEC(x.val_.val_, y.val_.val_);
  VEC g;
  a.d_.d_.grad(p, g);
  ASSERT_NEAR(third_deriv_aab, g[0], 0.03);
  ASSERT_NEAR(third_deriv_abb, g[1], 0.03);
}

struct log_falling_factorial_fun {
  template <typename T0, typename T1>
  inline typename boost::math::tools::promote_args<T0, T1>::type operator()(
      const T0 arg1, const T1 arg2) const {
    return log_falling_factorial(arg1, arg2);
  }
};

TEST(AgradFwdLogFallingFactorial, nan) {
  log_falling_factorial_fun log_falling_factorial_;
  test_nan_mix(log_falling_factorial_, 3.0, 5.0, false);
}

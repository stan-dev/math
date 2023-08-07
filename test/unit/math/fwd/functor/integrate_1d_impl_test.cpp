#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
/*
TEST(FwdFunctor, integrate_1d_t1) {
  using stan::math::fvar;
  using stan::math::integrate_1d_impl;
  using stan::math::EPSILON;

  double a = 0.2;
  double b = 0.7;
  fvar<double> theta = 0.75;
  theta.d_ = 1.0;
  const double relative_tolerance = std::sqrt(EPSILON);
  std::ostringstream *msgs = nullptr;

  auto func = [](const auto& x, const auto& xc,
                  std::ostream* msgs, const auto& theta) {
                    return stan::math::exp(x) + theta;
                  };

  fvar<double> res;
  res = integrate_1d_impl(func, a, b, relative_tolerance, msgs, theta);
  EXPECT_FLOAT_EQ(res.val(), 0.7923499493102901 + 0.5 * 0.75);
  EXPECT_FLOAT_EQ(res.d(), 0.5);

  res = integrate_1d_impl(func, 0.0, 0.0, relative_tolerance, msgs, theta);
  EXPECT_FLOAT_EQ(res.val(), 0);
  EXPECT_FLOAT_EQ(res.d(), 0);
}

TEST(FwdFunctor, integrate_1d_t2) {
  using stan::math::fvar;
  using stan::math::integrate_1d_impl;
  using stan::math::EPSILON;

  double a_dbl = 0.0;
  double b_dbl = 1.0;
  double theta_dbl = 0.5;
  fvar<double> a_fv(a_dbl, 1.0);
  fvar<double> b_fv(b_dbl, 1.0);
  fvar<double> theta_fv(theta_dbl, 1.0);
  const double relative_tolerance = std::sqrt(EPSILON);
  std::ostringstream *msgs = nullptr;

  auto func = [](const auto& x, const auto& xc, std::ostream* msgs,
                  const auto& theta) {
    return stan::math::exp(theta * stan::math::cos(2 * 3.141593 * x)) + theta;
  };

  fvar<double> res = integrate_1d_impl(func, a_fv, b_fv, relative_tolerance,
                                        msgs, theta_fv);
  EXPECT_FLOAT_EQ(res.val(), 1.56348343527304);
  EXPECT_FLOAT_EQ(res.d(),
                  1.25789445875152 -2.148721270700128 + 2.14872127069993);

  res = integrate_1d_impl(func, a_dbl, b_fv, relative_tolerance, msgs,
                          theta_fv);
  EXPECT_FLOAT_EQ(res.val(), 1.56348343527304);
  EXPECT_FLOAT_EQ(res.d(), 1.25789445875152 + 2.14872127069993);

  res = integrate_1d_impl(func, a_fv, b_dbl, relative_tolerance, msgs,
                          theta_fv);
  EXPECT_FLOAT_EQ(res.val(), 1.56348343527304);
  EXPECT_FLOAT_EQ(res.d(), 1.25789445875152 -2.148721270700128);

  res = integrate_1d_impl(func, a_fv, b_fv, relative_tolerance, msgs,
                          theta_dbl);
  EXPECT_FLOAT_EQ(res.val(), 1.56348343527304);
  EXPECT_FLOAT_EQ(res.d(), -2.148721270700128 + 2.14872127069993);
}
*/
TEST(FwdFunctor, integrate_1d_t3) {
  using stan::math::fvar;
  using stan::math::integrate_1d_impl;
  using stan::math::EPSILON;

  double a_dbl = 0.0;
  fvar<fvar<double>> a = a_dbl;
  a.val_.d_ = 1.0;
  a.d_.d_ = 1.0;

  double b_dbl = 1.0;

  fvar<fvar<double>> theta = 0.5;
  theta.d_ = 1.0;
  theta.d_.d_ = 1.0;

  const double relative_tolerance = std::sqrt(EPSILON);
  std::ostringstream *msgs = nullptr;

  auto func = [](const auto& x, const auto& xc, std::ostream* msgs,
                  const auto& theta) {
    return stan::math::exp(theta * stan::math::cos(2 * 3.141593 * x)) + theta;
  };

  fvar<fvar<double>> res = integrate_1d_impl(func, a_dbl, b_dbl,
                                              relative_tolerance, msgs, theta);
  EXPECT_FLOAT_EQ(res.val().val(), 1.56348343527304);
  EXPECT_FLOAT_EQ(res.val().d(), 1.25789445875152);

  // Second derivative not yet available
  EXPECT_FLOAT_EQ(res.d().val(), 1.25789445875152);
  EXPECT_FLOAT_EQ(res.d().d(), 0.5476947599595309);
}

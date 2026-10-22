#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>

TEST(mathMixScalFun, gammaQ) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::gamma_q(x1, x2);
  };
  stan::test::expect_ad(f, 0.5001, 1.0001);
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 1.0, 1.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  // this still fails forward mode; left regression test in rev/fun
  // stan::test::expect_value(f, 8.01006, 2.47579e+215);
}

/**
 * Forward-mode derivative with respect to the first argument at large z.
 * References are d/da Q(a, z) from mpmath at 80 digits.
 */
TEST(mathMixScalFun, gammaQ_fwd_shape_derivative_large_z) {
  using stan::math::fvar;
  struct Case {
    double a;
    double z;
    double dqda;
  };
  const Case cases[] = {
      {1.5, 300.0, 5.71509832058708679e-129},
      {3.0, 650.0, 6.01811954875449219e-277},
  };
  for (const auto& c : cases) {
    fvar<double> a = c.a;
    a.d_ = 1.0;
    fvar<double> z = c.z;
    const fvar<double> r = stan::math::gamma_q(a, z);
    ASSERT_TRUE(std::isfinite(r.d_)) << "a=" << c.a << " z=" << c.z;
    EXPECT_NEAR(c.dqda, r.d_, 1e-10 * std::fabs(c.dqda))
        << "a=" << c.a << " z=" << c.z;
  }
}

/**
 * Forward-mode derivative with respect to the second argument, at an a
 * above the range where tgamma(a) is finite. Reference from mpmath, closed
 * form -z^(a-1) e^(-z) / Gamma(a).
 */
TEST(mathMixScalFun, gammaQ_fwd_z_derivative_large_a) {
  using stan::math::fvar;
  fvar<double> a = 200.0;
  fvar<double> z = 150.0;
  z.d_ = 1.0;
  const fvar<double> r = stan::math::gamma_q(a, z);
  const double expected = -2.00507038377120400e-05;
  ASSERT_TRUE(std::isfinite(r.d_));
  EXPECT_NEAR(expected, r.d_, 1e-12 * std::fabs(expected));
}

TEST(MathFunctions, gammaQ_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::gamma_q;
    return gamma_q(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

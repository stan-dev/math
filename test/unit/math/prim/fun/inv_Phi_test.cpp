#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, inv_Phi) {
  using stan::math::inv_Phi;
  using stan::math::Phi;
  EXPECT_FLOAT_EQ(0.0, inv_Phi(0.5));
  double p = 0.123456789;
  EXPECT_FLOAT_EQ(p, Phi(inv_Phi(p)));
  p = 8e-311;
  EXPECT_FLOAT_EQ(p, Phi(inv_Phi(p)));
  p = 0.99;
  EXPECT_FLOAT_EQ(p, Phi(inv_Phi(p)));

  // breakpoints
  p = 0.02425;
  EXPECT_FLOAT_EQ(p, Phi(inv_Phi(p)));
  p = 0.97575;
  EXPECT_FLOAT_EQ(p, Phi(inv_Phi(p)));
}

TEST(MathFunctions, Equal){
  using stan::math::inv_Phi;
  // test output generated with R using qnorm
  EXPECT_NEAR(-2.247626755795137, inv_Phi(0.0123), 1e-15);
  EXPECT_NEAR(0.0582719627987081, inv_Phi(0.523234), 1e-15);
  EXPECT_NEAR(1.4284791211149008, inv_Phi(0.923423), 1e-15);
}

TEST(MathFunctions, inv_Phi_inf) {
  using stan::math::inv_Phi;
  double p = 7e-311;
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(inv_Phi(p), -inf);
  p = 1.0;
  EXPECT_EQ(inv_Phi(p), inf);
}
TEST(MathFunctions, inv_Phi_nan) {
  using stan::math::inv_Phi;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(inv_Phi(nan), std::domain_error);
  EXPECT_THROW(inv_Phi(-2.0), std::domain_error);
  EXPECT_THROW(inv_Phi(2.0), std::domain_error);
}

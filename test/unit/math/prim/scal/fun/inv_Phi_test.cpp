#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

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
TEST(MathFunctions, inv_Phi_inf) {
  using stan::math::inv_Phi;
  double p = 7e-311;
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(inv_Phi(p),-inf);
  p = 1.0;
  EXPECT_EQ(inv_Phi(p),inf);
}
TEST(MathFunctions, inv_Phi_nan) {
  using stan::math::inv_Phi;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(inv_Phi(nan), std::domain_error);
  EXPECT_THROW(inv_Phi(-2.0), std::domain_error);
  EXPECT_THROW(inv_Phi(2.0), std::domain_error);
}

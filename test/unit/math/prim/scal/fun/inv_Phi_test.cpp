#include <boost/math/special_functions/fpclassify.hpp>
#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>
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
}
TEST(MathFunctions, Phi_inf) {
  using stan::math::inv_Phi;
  double p = 7e-311;
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(inv_Phi(p),-inf);
  p = 0.99999999999999999999999999999999999999;
  EXPECT_EQ(inv_Phi(p),inf);
}
TEST(MathFunctions, Phi_nan) {
  using stan::math::inv_Phi;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(inv_Phi(nan), std::domain_error);
  EXPECT_THROW(inv_Phi(-2.0), std::domain_error);
  EXPECT_THROW(inv_Phi(2.0), std::domain_error);
}

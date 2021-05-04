#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(passesOnConvergentArgs, Check3F2Converges) {
  using stan::math::check_3F2_converges;
  const char* function = "check_3F2_converges";
  double a1 = 1.0;
  double a2 = 1.0;
  double a3 = 1.0;
  double b1 = 5.0;
  double b2 = 5.0;
  double z = 0.3;

  // in radius of convergence for z, other args don't matter
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 5.0;
  b2 = 5.0;
  // still in radius of convergence, ok
  z = 1.0;
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  // now in radius of convergences, but b1 is too small.
  b1 = 1.1;
  b2 = 1.1;
  z = 1.0;
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);

  // a1 is too big
  a1 = 40.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 10.0;
  b2 = 10.0;
  z = 1.0;
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);

  a1 = 5.0;
  a2 = 0.0;
  a3 = 1.0;
  b1 = 10.0;
  b2 = 10.0;
  z = 1.0;
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 5.0;
  b2 = 5.0;
  // outside of radius of convergence for current implementation.
  z = 1.3;
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 1.0;
  b2 = 1.0;
  // b1 is small, but z < 1 so we're ok.
  z = 0.99999999999;
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 1.0;
  b2 = 1.0;
  // checking negative z, this is fine
  z = -0.999999999999;
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 10.0;
  b2 = 10.0;
  // limits of range?
  z = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  b1 = 10.0;
  b2 = 10.0;
  // limits of range?
  z = -1.0 * std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);
  EXPECT_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z),
               std::domain_error);

  a1 = 1.0;
  a2 = 1.0;
  a3 = 1.0;
  // should be ok, underflow to zero (?)
  b1 = std::numeric_limits<double>::infinity();
  // should be ok, underflow to zero (?)
  b2 = std::numeric_limits<double>::infinity();
  z = 0.5;
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));
  EXPECT_NO_THROW(check_3F2_converges(function, a1, a2, a3, b1, b2, z));
}

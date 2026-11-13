#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log_bessel_i0e) {
  using stan::math::log_bessel_i0e;
  // z and log(I0(z)) - z, computed with mpmath (60 digits); the points include
  // both switches of Boost's approximation, at 7.75 and 500
  const double ref[][2]
      = {{1e-8, -9.999999975e-9},     {0.5, -0.4384502808145187},
         {7.7, -1.9220767903844019},  {7.75, -1.9254352770188191},
         {100, -3.2202673100574163},  {499.5, -4.0254918309739383},
         {500, -4.0259923318933035},  {1e4, -5.524096218567699},
         {1e12, -14.734449091168822}, {1e300, -346.30670248231153}};
  for (const auto& r : ref) {
    EXPECT_NEAR(r[1], log_bessel_i0e(r[0]), 1e-15 * std::fabs(r[1]));
  }
  EXPECT_EQ(0.0, log_bessel_i0e(0.0));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            log_bessel_i0e(std::numeric_limits<double>::infinity()));
  // integer argument
  EXPECT_NEAR(-3.2202673100574163, log_bessel_i0e(100), 1e-15 * 3.3);

  Eigen::VectorXd z(2);
  z << 0.5, 1e12;
  Eigen::VectorXd v = log_bessel_i0e(z);
  EXPECT_NEAR(-0.4384502808145187, v(0), 1e-15);
  EXPECT_NEAR(-14.734449091168822, v(1), 1e-14);

  EXPECT_THROW(log_bessel_i0e(-1.0), std::domain_error);
  EXPECT_THROW(log_bessel_i0e(std::numeric_limits<double>::quiet_NaN()),
               std::domain_error);
}

#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, lambert_w1) {
  using stan::math::lambert_w1;
  using stan::math::lambert_w0;
  printf("%f %s", lambert_w1(-0.2), "\n");
  printf("%f %s", lambert_w1(-0.3), "\n");

  printf("%f %s", lambert_w0(1), "\n");
  printf("%f %s", lambert_w0(0), "\n");
  printf("%f %s", lambert_w0(-1), "\n");

}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

// This does not match the example at
// https://en.cppreference.com/w/cpp/numeric/complex/proj
TEST(AgradRev, proj_test) {
  std::complex<stan::math::var> z1(1, 2);
  auto f = proj(z1);
  EXPECT_TRUE(f == z1) << f << std::endl;

  using stan::math::positive_infinity;
  std::complex<stan::math::var> z2(positive_infinity(), -1);
  f = proj(z2);
  EXPECT_TRUE(f == std::complex<stan::math::var>(positive_infinity(), -0))
      << f << std::endl;

  using stan::math::negative_infinity;
  std::complex<stan::math::var> z3(0, negative_infinity());
  f = proj(z3);
  EXPECT_TRUE(f == std::complex<stan::math::var>(positive_infinity(), -0))
      << f << std::endl;
}

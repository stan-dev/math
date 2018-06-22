#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, polar_test) {
  std::complex<stan::math::var> z = std::polar(1, 0);
  auto f = arg(z);
  EXPECT_EQ(f.val(), 0.0);
}

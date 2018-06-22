#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, arg_test) {
  std::complex<stan::math::var> z
      = std::complex<stan::math::var>(1, stan::math::pi());
  auto f = arg(z);
  using stan::math::atan2;
  EXPECT_EQ(f.val(), atan2(std::imag(z), std::real(z)));
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, conj_test) {
  std::complex<stan::math::var> z(1, 2);
  auto f = conj(z);
  EXPECT_EQ(real(f).val(), 1);
  EXPECT_EQ(imag(f).val(), -2);
}

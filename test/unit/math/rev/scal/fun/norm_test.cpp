#include <stan/math/rev/scal.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, norm_test) {
  std::complex<stan::math::var> z = std::complex<stan::math::var>(3, 4);
  auto f = norm(z);
  EXPECT_EQ(f.val(), 25);

  AVEC x = createAVEC(imag(z));
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(g[0], 8);
}

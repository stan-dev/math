#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, asin_var) {
  AVAR a = 0.68;
  AVAR f = asin(a);
  EXPECT_FLOAT_EQ(asin(0.68), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / sqrt(1.0 - (0.68 * 0.68)), g[0]);
}

TEST(AgradRev, asin_1) {
  AVAR a = 1;
  AVAR f = asin(a);
  EXPECT_FLOAT_EQ((1.57079632679), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / sqrt(1.0 - (1 * 1)), g[0]);
}

TEST(AgradRev, asin_neg_1) {
  AVAR a = -1;
  AVAR f = asin(a);
  EXPECT_FLOAT_EQ((-1.57079632679), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / sqrt(1.0 - (-1 * -1)), g[0]);
}

TEST(AgradRev, asin_out_of_bounds1) {
  AVAR a = 1.0 + stan::math::EPSILON;
  AVAR f = asin(a);
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_TRUE(std::isnan(asin(a)));
  EXPECT_EQ(g.size(), 1);
  EXPECT_TRUE(std::isnan(g[0]));
}
TEST(AgradRev, asin_out_of_bounds2) {
  AVAR a = -1.0 - stan::math::EPSILON;
  AVAR f = asin(a);
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_TRUE(std::isnan(asin(a)));
  EXPECT_EQ(g.size(), 1);
  EXPECT_TRUE(std::isnan(g[0]));
}

struct asin_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return asin(arg1);
  }
};

TEST(AgradRev, asin_NaN) {
  asin_fun asin_;
  test_nan(asin_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.68;
  test::check_varis_on_stack(stan::math::asin(a));
}

TEST(AgradRev, complex) {
  stan::math::var x = 2.0 / 3;

  double h = 1e-8;
  std::complex<stan::math::var> z(x, h);
  auto f = asin(z);

  AVEC v = createAVEC(real(z));
  VEC g;
  real(f).grad(v, g);
  EXPECT_FLOAT_EQ(g[0], imag(f).val() / h);
}

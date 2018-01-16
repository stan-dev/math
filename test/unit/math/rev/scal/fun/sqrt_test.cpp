#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, sqrt_a) {
  AVAR a(5.0);
  AVAR f = sqrt(a);
  EXPECT_FLOAT_EQ(sqrt(5.0), f.val());
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ((1.0 / 2.0) * pow(5.0, -0.5), g[0]);
}

TEST(AgradRev, sqrt_neg) {
  AVAR a = 0.0 - stan::math::EPSILON;
  EXPECT_TRUE(std::isnan(sqrt(a)));

  a = -100;
  EXPECT_TRUE(std::isnan(sqrt(a)));
}

TEST(AgradRev, sqrt_inf) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = inf;
  AVAR f = sqrt(a);
  EXPECT_FLOAT_EQ(inf, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.0, g[0]);
}

TEST(AgradRev, sqrt_zero) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a(0.0);
  AVAR f = sqrt(a);
  EXPECT_FLOAT_EQ(0.0, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(inf, g[0]);
}

struct sqrt_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return sqrt(arg1);
  }
};

TEST(AgradRev, sqrt_NaN) {
  sqrt_fun sqrt_;
  test_nan(sqrt_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a(5.0);
  test::check_varis_on_stack(stan::math::sqrt(a));
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, atan_1) {
  AVAR a = 1;
  AVAR f = atan(a);
  EXPECT_FLOAT_EQ((.78539816339), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / (1.0 + (1 * 1)), g[0]);
}

TEST(AgradRev, atan_neg_1) {
  AVAR a = -1;
  AVAR f = atan(a);
  EXPECT_FLOAT_EQ((-.78539816339), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / (1.0 + (-1 * -1)), g[0]);
}

TEST(AgradRev, atan_boundry) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = inf;
  AVAR f = atan(a);
  EXPECT_FLOAT_EQ(1.5707964, f.val());
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.0, g[0]);

  AVAR b = -inf;
  AVAR e = atan(b);
  EXPECT_FLOAT_EQ(-1.5707964, e.val());
  AVEC y = createAVEC(b);
  VEC h;
  e.grad(y, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
}

struct atan_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return atan(arg1);
  }
};

TEST(AgradRev, atan_NaN) {
  atan_fun atan_;
  test_nan(atan_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1;
  test::check_varis_on_stack(stan::math::atan(a));
}

TEST(AgradRev, complex) {
  stan::math::var x = 2.0 / 3;

  double h = 1e-8;
  std::complex<stan::math::var> z(x, h);
  auto f = atan(z);

  AVEC v = createAVEC(real(z));
  VEC g;
  real(f).grad(v, g);
  EXPECT_FLOAT_EQ(g[0], imag(f).val() / h);
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, abs_var) {
  AVAR a = 0.68;
  AVAR f = abs(a);
  EXPECT_FLOAT_EQ(0.68, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0, g[0]);
}

TEST(AgradRev, abs_var_2) {
  AVAR a = -0.68;
  AVAR f = abs(a);
  EXPECT_FLOAT_EQ(0.68, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(-1.0, g[0]);
}

TEST(AgradRev, abs_var_3) {
  AVAR a = 0.0;
  AVAR f = abs(a);
  EXPECT_FLOAT_EQ(0.0, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_EQ(1U, g.size());
  EXPECT_FLOAT_EQ(0.0, g[0]);
}

TEST(AgradRev, abs_inf) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = inf;
  AVAR f = abs(a);
  EXPECT_FLOAT_EQ(inf, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0, g[0]);
}

TEST(AgradRev, abs_neg_inf) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = -inf;
  AVAR f = abs(a);
  EXPECT_FLOAT_EQ(inf, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(-1.0, g[0]);
}

struct abs_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return abs(arg1);
  }
};

TEST(AgradRev, abs_NaN) {
  abs_fun abs_;
  test_nan(abs_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.68;
  test::check_varis_on_stack(stan::math::abs(a));
}

TEST(AgradRev, abs_complex) {
  std::complex<stan::math::var> z = std::complex<stan::math::var>(3, 4);
  auto f = abs(z);
  EXPECT_EQ(5, f.val());
  AVEC x = createAVEC(real(z));
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.6, g[0]);
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, pow_var_var) {
  AVAR a(3.0);
  AVAR b(4.0);
  AVAR f = pow(a, b);
  EXPECT_FLOAT_EQ(81.0, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(4.0 * pow(3.0, 4.0 - 1.0), g[0]);
  EXPECT_FLOAT_EQ(log(3.0) * pow(3.0, 4.0), g[1]);
}

TEST(AgradRev, pow_var_double) {
  AVAR a(4.0);
  double b = 4.0;
  AVAR f = pow(a, b);
  EXPECT_FLOAT_EQ(256.0, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(4.0 * pow(4.0, 4.0 - 1.0), g[0]);

  b = 2.0;
  f = pow(a, b);
  EXPECT_FLOAT_EQ(16.0, f.val());

  b = 0.5;
  f = pow(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  b = 1.0;
  f = pow(a, b);
  EXPECT_FLOAT_EQ(a.val(), f.val());

  b = -0.5;
  f = pow(a, b);
  EXPECT_FLOAT_EQ(0.5, f.val());

  b = -2.0;
  f = pow(a, b);
  EXPECT_FLOAT_EQ(1 / 16.0, f.val());
}

TEST(AgradRev, pow_double_var) {
  double a = 3.0;
  AVAR b(4.0);
  AVAR f = pow(a, b);
  EXPECT_FLOAT_EQ(81.0, f.val());

  AVEC x = createAVEC(b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(log(3.0) * pow(3.0, 4.0), g[0]);
}

TEST(AgradRev, pow_boundry) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = inf;
  AVAR b = 5;
  AVAR f = pow(a, b);
  EXPECT_FLOAT_EQ(inf, f.val());
  AVAR g = pow(b, a);
  EXPECT_FLOAT_EQ(inf, g.val());

  AVAR c = -inf;
  AVAR d = 6;
  AVAR h = pow(c, b);
  EXPECT_FLOAT_EQ(-inf, h.val());
  AVAR i = pow(c, d);
  EXPECT_FLOAT_EQ(inf, i.val());

  AVAR j = pow(b, c);
  EXPECT_FLOAT_EQ(0.0, j.val());
}

struct pow_fun {
  template <typename T0, typename T1>
  inline typename stan::return_type<T0, T1>::type operator()(
      const T0& arg1, const T1& arg2) const {
    return pow(arg1, arg2);
  }
};

TEST(AgradRev, pow_nan) {
  pow_fun pow_;
  test_nan(pow_, 3.0, 5.0, false, true);
  test_nan(pow_, 0.0, 5.0, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a(3.0);
  AVAR b(4.0);
  test::check_varis_on_stack(stan::math::pow(a, b));
  test::check_varis_on_stack(stan::math::pow(3.0, b));
  test::check_varis_on_stack(stan::math::pow(a, 4.0));
}

// this test used to not build with g++-4.9 and lower
TEST(AgradRev, complex) {
  std::complex<stan::math::var> i(0, 1);
  auto f = pow(i, i);
  EXPECT_EQ(real(f).val(), exp(-stan::math::pi() / 2));
}

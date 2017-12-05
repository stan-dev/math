#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, fmin_vv) {
  AVAR a = 1.3;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(1.3, f.val());

  AVEC x = createAVEC(a, b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(1.0, grad_f[0]);
  EXPECT_FLOAT_EQ(0.0, grad_f[1]);
}

TEST(AgradRev, fmin_vv_2) {
  AVAR a = 2.3;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(a, b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
  EXPECT_FLOAT_EQ(1.0, grad_f[1]);
}

TEST(AgradRev, fmin_vv_3) {
  AVAR a = 2.0;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(a, b);
  VEC grad_f;
  f.grad(x, grad_f);
  // arbitrary, but documented this way
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
  EXPECT_FLOAT_EQ(1.0, grad_f[1]);
}

TEST(AgradRev, fmin_vd) {
  AVAR a = 1.3;
  double b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(1.3, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(1.0, grad_f[0]);
}

TEST(AgradRev, fmin_vd_2) {
  AVAR a = 2.3;
  double b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

TEST(AgradRev, fmin_vd_3) {
  AVAR a = 2.0;
  double b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(1.0, grad_f[0]);
}

TEST(AgradRev, fmin_dv) {
  double a = 1.3;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(1.3, f.val());

  AVEC x = createAVEC(b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

TEST(AgradRev, fmin_dv_2) {
  double a = 2.3;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(1.0, grad_f[0]);
}

TEST(AgradRev, fmin_dv_3) {
  double a = 2.0;
  AVAR b = 2.0;
  AVAR f = fmin(a, b);
  EXPECT_FLOAT_EQ(2.0, f.val());

  AVEC x = createAVEC(b);
  VEC grad_f;
  f.grad(x, grad_f);
  // arbitrary, but doc this way
  EXPECT_FLOAT_EQ(1.0, grad_f[0]);
}

struct fmin_fun {
  template <typename T0, typename T1>
  inline
  typename stan::return_type<T0, T1>::type
  operator()(const T0& arg1,
             const T1& arg2) const {
    return fmin(arg1, arg2);
  }
};

TEST(AgradRev, fmin_nan) {
  fmin_fun fmin_;
  double nan = std::numeric_limits<double>::quiet_NaN();
  test_nan_vv(fmin_, nan, nan, false, true);
  test_nan_dv(fmin_, nan, nan, false, true);
  test_nan_vd(fmin_, nan, nan, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  AVAR b = 2.0;
  test::check_varis_on_stack(stan::math::fmin(a, b));
  test::check_varis_on_stack(stan::math::fmin(a, 2.0));
  test::check_varis_on_stack(stan::math::fmin(1.3, b));
}

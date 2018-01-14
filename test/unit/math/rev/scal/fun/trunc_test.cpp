#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, trunc) {
  AVAR a = 1.2;
  AVAR f = stan::math::trunc(a);
  EXPECT_FLOAT_EQ(1.0, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

TEST(AgradRev, trunc_2) {
  AVAR a = -1.2;
  AVAR f = stan::math::trunc(a);
  EXPECT_FLOAT_EQ(-1.0, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

struct trunc_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return trunc(arg1);
  }
};

TEST(AgradRev, trunc_NaN) {
  trunc_fun trunc_;
  test_nan(trunc_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = -1.2;
  test::check_varis_on_stack(stan::math::trunc(a));
}

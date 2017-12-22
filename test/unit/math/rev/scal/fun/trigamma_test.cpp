#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, trigamma) {
  AVAR a = 0.5;
  AVAR f = stan::math::trigamma(a);
  EXPECT_FLOAT_EQ(4.9348022005446793094, f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(-16.8288, grad_f[0]);
}

struct trigamma_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return stan::math::trigamma(arg1);
  }
};

TEST(AgradRev, trigamma_NaN) {
  trigamma_fun trigamma_;
  test_nan(trigamma_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.5;
  test::check_varis_on_stack(stan::math::trigamma(a));
}

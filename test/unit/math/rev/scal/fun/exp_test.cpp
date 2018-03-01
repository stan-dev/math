#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, exp_a) {
  AVAR a(6.0);
  // mix exp() functs w/o namespace
  AVAR f = exp(a);
  EXPECT_FLOAT_EQ(exp(6.0), f.val());
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(exp(6.0), g[0]);
}

struct exp_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return exp(arg1);
  }
};

TEST(AgradRev, exp_NaN) {
  exp_fun exp_;
  test_nan(exp_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a(6.0);
  test::check_varis_on_stack(stan::math::exp(a));
}

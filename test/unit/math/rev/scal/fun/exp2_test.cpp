#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, exp2) {
  AVAR a = 1.3;
  AVAR f = stan::math::exp2(a);
  EXPECT_FLOAT_EQ(std::pow(2.0, 1.3), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(std::pow(2.0, 1.3) * std::log(2.0), grad_f[0]);

  a = std::numeric_limits<AVAR>::infinity();
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),
                  stan::math::exp2(a).val());
}

struct exp2_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return exp2(arg1);
  }
};

TEST(AgradRev, exp2_NaN) {
  exp2_fun exp2_;
  test_nan(exp2_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  test::check_varis_on_stack(stan::math::exp2(a));
}

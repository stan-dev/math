#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, log2) {
  AVAR a = 3.0;
  AVAR f = stan::math::log2(a);
  EXPECT_FLOAT_EQ(std::log(3.0)/std::log(2.0), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(1.0 / 3.0 / std::log(2.0), grad_f[0]);

  a = std::numeric_limits<AVAR>::infinity();
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),
                  stan::math::log2(a).val());
}

struct log2_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log2(arg1);
  }
};

TEST(AgradRev, log2_NaN) {
  log2_fun log2_;
  test_nan(log2_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 3.0;
  test::check_varis_on_stack(stan::math::log2(a));
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

void test_logit(double u) {
  using stan::math::log;
  using stan::math::logit;
  AVAR u1 = u;
  AVEC uv1 = createAVEC(u1);
  AVAR f1 = log(u1 / (1 - u1));
  VEC grad_f1;
  f1.grad(uv1, grad_f1);
  double g1 = grad_f1[0];

  AVAR u2 = u;
  AVEC uv2 = createAVEC(u2);
  AVAR f2 = logit(u2);
  VEC grad_f2;
  f2.grad(uv2, grad_f2);
  double g2 = grad_f2[0];

  EXPECT_FLOAT_EQ(log(u / (1 - u)), f2.val());
  EXPECT_FLOAT_EQ(g1, g2);
}

TEST(AgradRev, logitDeriv) {
  test_logit(0.3);
  test_logit(0.5);
}

struct logit_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return stan::math::logit(arg1);
  }
};

TEST(AgradRev, inv_logit_NaN) {
  logit_fun logit_;
  test_nan(logit_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR u = 0.3;
  test::check_varis_on_stack(stan::math::logit(u));
}

#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, tgamma) {
  AVAR a = 3.5;
  AVAR f = stan::math::tgamma(a);
  EXPECT_FLOAT_EQ(stan::math::tgamma(3.5), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(boost::math::digamma(3.5) * stan::math::tgamma(3.5),
                  grad_f[0]);
}

struct tgamma_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return tgamma(arg1);
  }
};

TEST(AgradRev, tgamma_NaN) {
  tgamma_fun tgamma_;
  test_nan(tgamma_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 3.5;
  test::check_varis_on_stack(stan::math::tgamma(a));
}

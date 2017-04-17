#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <stdexcept>

TEST(AgradRev,log1p) {
  AVAR a = 0.1;
  AVAR f = stan::math::log1p(a);
  EXPECT_FLOAT_EQ(log(1 + 0.1), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x,grad_f);
  EXPECT_FLOAT_EQ(1.0 / (1.0 + 0.1), grad_f[0]);
}

TEST(AgradRevLog1p, excepts) {
  using stan::math::log1p;
  AVAR a = -2;
  EXPECT_THROW(log1p(a), std::domain_error);
}

TEST(AgradRevLog1p, overflows) {
  using stan::math::log1p;
  AVAR b = -1;
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::log1p(b).val());
}

struct log1p_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log1p(arg1);
  }
};

TEST(AgradRev,log1p_NaN) {
  log1p_fun log1p_;
  test_nan(log1p_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.1;
  test::check_varis_on_stack(stan::math::log1p(a));
}

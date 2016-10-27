#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,log1m) {
  using stan::math::log1m;
  AVAR a = 0.1;
  AVAR f = log1m(a);
  EXPECT_FLOAT_EQ(log(1 - 0.1), f.val());
  
  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x,grad_f);
  EXPECT_FLOAT_EQ(-1 / (1 - 0.1), grad_f[0]);
}
TEST(AgradRev,excepts) {
  using stan::math::log1m;
  EXPECT_THROW(log1m(AVAR(10)), std::domain_error);
}

TEST(MathFunctions, log1m_inf_return) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::log1m(1));
}

struct log1m_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log1m(arg1);
  }
};

TEST(AgradRev,log1m_NaN) {
  log1m_fun log1m_;
  test_nan(log1m_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.1;
  test::check_varis_on_stack(stan::math::log1m(a));
}

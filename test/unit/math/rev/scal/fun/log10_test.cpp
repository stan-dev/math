#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,log10_a) {
  AVAR a(5.0);
  AVAR f = log10(a); 
  EXPECT_FLOAT_EQ(log10(5.0),f.val());
  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(1.0/(log(10.0) * 5.0),g[0]);
}

struct log10_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log10(arg1);
  }
};

TEST(AgradRev,log10_NaN) {
  log10_fun log10_;
  test_nan(log10_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  stan::math::var x = 4.0;
  test::check_varis_on_stack(stan::math::log10(x));
}

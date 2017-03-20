#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/zeta.hpp>

TEST(AgradRev,digamma) {
  AVAR a = 0.5;
  AVAR f = digamma(a);
  EXPECT_FLOAT_EQ(boost::math::digamma(0.5),f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x,grad_f);
  EXPECT_FLOAT_EQ(4.9348022005446793094, grad_f[0]);
}  

struct digamma_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return digamma(arg1);
  }
};

TEST(AgradRev,digamma_NaN) {
  digamma_fun digamma_;
  test_nan(digamma_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.5;
  test::check_varis_on_stack(stan::math::digamma(a));
}

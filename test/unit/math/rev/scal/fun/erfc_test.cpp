#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,erfc) {
  AVAR a = 1.3;
  AVAR f = erfc(a);
  EXPECT_FLOAT_EQ(::erfc(1.3), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x,grad_f);
  EXPECT_FLOAT_EQ(-2.0 / std::sqrt(boost::math::constants::pi<double>()) * std::exp(- 1.3 * 1.3), grad_f[0]);
}

struct erfc_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return erfc(arg1);
  }
};

TEST(AgradRev,erfc_NaN) {
  erfc_fun erfc_;
  test_nan(erfc_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  test::check_varis_on_stack(stan::math::erfc(a));
}

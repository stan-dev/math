#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,modified_bessel_second_kind_int_var) {
  int a(1);
  AVAR b(4.0);
  AVAR f = stan::math::modified_bessel_second_kind(a,b);
  EXPECT_FLOAT_EQ(0.01248349888726843147038417998080606848384,f.val());

  AVEC x = createAVEC(a,b);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(0,g[0]);
  EXPECT_FLOAT_EQ(-0.01428055080767013213734124, g[1]);

  a = 1;
  b = -3.0;
  EXPECT_THROW(stan::math::modified_bessel_second_kind(a,b), std::domain_error);

  a = -1;
  EXPECT_THROW(stan::math::modified_bessel_second_kind(a,b), std::domain_error);
}

struct modified_bessel_second_kind_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return modified_bessel_second_kind(2,arg1);
  }
};

TEST(AgradRev,modified_bessel_second_kind_NaN) {
  modified_bessel_second_kind_fun modified_bessel_second_kind_;
  test_nan(modified_bessel_second_kind_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  int a(1);
  AVAR b(4.0);
  test::check_varis_on_stack(stan::math::modified_bessel_second_kind(a, b));
}

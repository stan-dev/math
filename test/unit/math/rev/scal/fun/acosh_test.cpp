#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <cmath>
#include <limits>

TEST(AgradRev,acosh_val) {
  using stan::math::acosh;
  using std::sqrt;
  AVAR a = 1.3;
  AVAR f = acosh(a);
  EXPECT_FLOAT_EQ(acosh(1.3), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1 / sqrt(1.3 * 1.3  - 1.0), g[0]);
}

TEST(AgradRev,acosh_1) {
  using stan::math::acosh;
  using std::sqrt;
  AVAR a = 1.0;
  AVAR f = acosh(a);
  EXPECT_FLOAT_EQ(0.0, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1 / sqrt(-1 * -1 - 1), g[0]);
}

TEST(MathFunctions, acosh_exception) {
  using stan::math::var;
  using stan::math::acosh;
  EXPECT_THROW(acosh(var(0.5)), std::domain_error);
}

TEST(AgradRevAcosh, overflows) {
  using stan::math::acosh;
  AVAR b = std::numeric_limits<double>::infinity();
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::acosh(b).val());
}


struct acosh_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    using stan::math::acosh;
    return acosh(arg1);
  }
};

TEST(AgradRev,acosh_NaN) {
  acosh_fun acosh_;
  test_nan(acosh_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  test::check_varis_on_stack(stan::math::acosh(a));
}

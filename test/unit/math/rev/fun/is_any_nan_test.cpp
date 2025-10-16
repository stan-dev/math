#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(MathFunctions, is_any_nan_variadic_rev) {
  using stan::math::is_any_nan;

  double dbl_inf = std::numeric_limits<double>::infinity();
  double dbl_nan = std::numeric_limits<double>::quiet_NaN();

  stan::math::var var_nan(dbl_nan);
  stan::math::var var_inf(dbl_inf);

  stan::math::var a = 7.0;
  stan::math::var b = 2.0;

  EXPECT_TRUE(is_any_nan(var_nan, b, 6, 5.0));
  EXPECT_TRUE(is_any_nan(var_nan, var_inf, 6, 5.0));
  EXPECT_TRUE(is_any_nan(dbl_nan, a, b, 6, 5.0));

  EXPECT_FALSE(is_any_nan(dbl_inf, b, 6, var_inf));
  EXPECT_FALSE(is_any_nan(a, b, 6, 7.0));
}

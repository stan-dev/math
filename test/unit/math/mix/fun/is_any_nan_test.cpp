#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(AgradMixIsAnyNan, Fvar) {
  using stan::math::fvar;
  using stan::math::is_any_nan;
  using stan::math::var;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<var> fv_nan(nan, nan);
  fvar<var> fv_max(max, max);
  fvar<var> fv_min(min, min);
  fvar<var> fv_dbl(0.5, 1.0);
  fvar<var> fv_inf(infinity, infinity);

  EXPECT_TRUE(is_any_nan(fv_nan, fv_dbl, 7.0));
  EXPECT_TRUE(is_any_nan(nan, fv_dbl, fv_min));
  EXPECT_TRUE(is_any_nan(fv_dbl, 1.5, 6, nan));

  EXPECT_FALSE(is_any_nan(fv_inf, fv_dbl, 7.0, min));
  EXPECT_FALSE(is_any_nan(max, infinity, fv_min, 8));
  EXPECT_FALSE(is_any_nan(fv_dbl, fv_inf, 6, min));
}

TEST(AgradMixIsAnyNan, FvarFvar) {
  using stan::math::fvar;
  using stan::math::is_any_nan;
  using stan::math::var;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<fvar<var>> ffv_nan(nan, nan);
  fvar<fvar<var>> ffv_max(max, max);
  fvar<fvar<var>> ffv_min(min, min);
  fvar<fvar<var>> ffv_dbl(0.5, 1.0);
  fvar<fvar<var>> ffv_inf(infinity, infinity);

  EXPECT_TRUE(is_any_nan(ffv_nan, ffv_dbl, 7.0));
  EXPECT_TRUE(is_any_nan(nan, ffv_dbl, ffv_min));
  EXPECT_TRUE(is_any_nan(ffv_dbl, 1.5, 6, nan));

  EXPECT_FALSE(is_any_nan(ffv_inf, ffv_dbl, 7.0, min));
  EXPECT_FALSE(is_any_nan(max, infinity, ffv_min, 8));
  EXPECT_FALSE(is_any_nan(ffv_dbl, ffv_inf, 6, min));
}

#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(AgradFwdIsAnyNan, Fvar) {
  using stan::math::fvar;
  using stan::math::is_any_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<double> fd_nan(nan, nan);
  fvar<double> fd_max(max, max);
  fvar<double> fd_min(min, min);
  fvar<double> fd_dbl(0.5, 1.0);
  fvar<double> fd_inf(infinity, infinity);

  EXPECT_TRUE(is_any_nan(fd_nan, fd_dbl, 7.0));
  EXPECT_TRUE(is_any_nan(nan, fd_dbl, fd_min));
  EXPECT_TRUE(is_any_nan(fd_dbl, 1.5, 6, nan));

  EXPECT_FALSE(is_any_nan(fd_inf, fd_dbl, 7.0, min));
  EXPECT_FALSE(is_any_nan(max, infinity, fd_min, 8));
  EXPECT_FALSE(is_any_nan(fd_dbl, fd_inf, 6, min));
}

TEST(AgradFwdIsAnyNan, FvarFvar) {
  using stan::math::fvar;
  using stan::math::is_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<fvar<double>> ffd_nan(nan, nan);
  fvar<fvar<double>> ffd_max(max, max);
  fvar<fvar<double>> ffd_min(min, min);
  fvar<fvar<double>> ffd_dbl(0.5, 1.0);
  fvar<fvar<double>> ffd_inf(infinity, infinity);

  EXPECT_TRUE(is_any_nan(ffd_nan, ffd_dbl, 7.0));
  EXPECT_TRUE(is_any_nan(nan, ffd_dbl, ffd_min));
  EXPECT_TRUE(is_any_nan(ffd_dbl, 1.5, 6, nan));

  EXPECT_FALSE(is_any_nan(ffd_inf, ffd_dbl, 7.0, min));
  EXPECT_FALSE(is_any_nan(max, infinity, ffd_min, 8));
  EXPECT_FALSE(is_any_nan(ffd_dbl, ffd_inf, 6, min));
}

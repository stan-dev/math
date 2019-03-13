#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandling, checkMatchingSizes) {
  std::vector<double> a;
  std::vector<double> b;

  EXPECT_NO_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b));

  a = {3, 3};
  EXPECT_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b),
      std::invalid_argument);

  b = {3, 3};
  EXPECT_NO_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b));
}

TEST(ErrorHandling, checkMatchingSizes_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> a;
  std::vector<double> b;
  EXPECT_NO_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b));

  a = {nan, nan};
  EXPECT_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b),
      std::invalid_argument);

  b = {nan, nan};
  EXPECT_NO_THROW(
      stan::math::check_matching_sizes("checkMatchingSizes", "a", a, "b", b));
}

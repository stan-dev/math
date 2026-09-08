#include <stan/math/prim/err/check_consistent_sizes.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(ErrorHandling, checkConsistentSizesUnnamedArguments) {
  using stan::math::check_consistent_sizes;

  std::vector<double> x1(3);
  std::vector<int> x2(3);
  std::vector<double> x3(2);

  EXPECT_NO_THROW(check_consistent_sizes("check_consistent_sizes", x1, x2));
  EXPECT_THROW_MSG_WITH_COUNT(
      check_consistent_sizes("check_consistent_sizes", x1, x2, x3),
      std::invalid_argument, "arg1", 1);
  EXPECT_THROW_MSG_WITH_COUNT(
      check_consistent_sizes("check_consistent_sizes", x1, x2, x3),
      std::invalid_argument, "arg3", 1);
  EXPECT_NO_THROW(check_consistent_sizes("check_consistent_sizes", x1));
}

TEST(ErrorHandling, checkConsistentSizesNamedArgumentsStillResolve) {
  using stan::math::check_consistent_sizes;

  std::vector<double> x1(3);
  std::vector<double> x2(2);

  EXPECT_THROW_MSG(
      check_consistent_sizes("check_consistent_sizes", "x1", x1, "x2", x2),
      std::invalid_argument, "x2");
}

#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandling, checkMatchingSizes) {
  std::vector<double> a;
  std::vector<double> b;

  EXPECT_NO_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                   "b", b));

  a.push_back(3.0);
  a.push_back(3.0);

  EXPECT_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                "b", b),
               std::invalid_argument);

  b.push_back(3.0);
  b.push_back(3.0);
  EXPECT_NO_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                   "b", b));
}

TEST(ErrorHandling, checkMatchingSizes_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> a;
  std::vector<double> b;
  EXPECT_NO_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                   "b", b));

  a.push_back(nan);
  a.push_back(nan);
  EXPECT_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                "b", b),
               std::invalid_argument);

  b.push_back(nan);
  b.push_back(nan);
  EXPECT_NO_THROW(stan::math::check_matching_sizes("checkMatchingSizes", "a", a,
                                                   "b", b));
}


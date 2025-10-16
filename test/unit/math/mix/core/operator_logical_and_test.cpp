#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

void test_logical_and(double x, double y) {
  using stan::math::fvar;
  using stan::math::var;

  EXPECT_EQ(x && y, fvar<double>(x) && fvar<double>(y));
  EXPECT_EQ(x && y, fvar<double>(x) && y);
  EXPECT_EQ(x && y, x && fvar<double>(y));

  EXPECT_EQ(x && y, fvar<fvar<double> >(x) && fvar<fvar<double> >(y));
  EXPECT_EQ(x && y, fvar<fvar<double> >(x) && y);
  EXPECT_EQ(x && y, x && fvar<fvar<double> >(y));

  EXPECT_EQ(x && y, fvar<var>(x) && fvar<var>(y));
  EXPECT_EQ(x && y, fvar<var>(x) && y);
  EXPECT_EQ(x && y, x && fvar<var>(y));

  EXPECT_EQ(x && y, fvar<fvar<var> >(x) && fvar<fvar<var> >(y));
  EXPECT_EQ(x && y, fvar<fvar<var> >(x) && y);
  EXPECT_EQ(x && y, x && fvar<fvar<var> >(y));
}

TEST(AgradMix, logical_and) {
  std::vector<double> xs;
  xs.push_back(6.1);
  xs.push_back(6.1);
  xs.push_back(0);
  xs.push_back(-13.2);
  xs.push_back(std::numeric_limits<double>::infinity());
  xs.push_back(-std::numeric_limits<double>::infinity());
  xs.push_back(std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < xs.size(); ++i)
    for (size_t j = 0; j < xs.size(); ++j)
      test_logical_and(xs[i], xs[j]);
}

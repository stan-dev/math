#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

template <typename T_N, typename T_n>
void test_choose(const T_N& N, const T_n& n) {
  using stan::math::choose;
  EXPECT_FLOAT_EQ(std::exp(lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1)),
                  choose(N, n));
}


TEST(MathFunctions, choose) {
  using stan::math::choose;
  EXPECT_EQ(1, choose(2, 2));
  EXPECT_EQ(2, choose(2, 1));
  EXPECT_EQ(3, choose(3, 1));
  EXPECT_EQ(3, choose(3, 2));
  EXPECT_EQ(0, choose(2, 3));
  for (int n = 0; n < 32; ++n) {
    test_choose(32, n);
  }
  EXPECT_THROW(choose(36, 18), std::domain_error);
}

TEST(MathFunctions, choose_nan) {
  using stan::math::choose;
  int nan = std::numeric_limits<int>::quiet_NaN() - 1;
  // quiet_NaN() returns 0 which would otherwise be valid
  EXPECT_THROW(choose(2, nan), std::domain_error);
  EXPECT_THROW(choose(nan, 2), std::domain_error);
  EXPECT_THROW(choose(nan, nan), std::domain_error);

}


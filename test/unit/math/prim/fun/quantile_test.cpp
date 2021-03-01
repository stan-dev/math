#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

template <typename T>
void test_quantile_double() {
  using stan::math::index_type_t;
  using stan::math::quantile;

  T c(1);
  c[0] = 1.7;
  EXPECT_EQ(quantile(c, 0), 1.7);
  EXPECT_EQ(quantile(c, 1), 1.7);
  EXPECT_EQ(quantile(c, 0.33), 1.7);
  EXPECT_EQ(quantile(c, 0.68), 1.7);

  T v(5);
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  v[3] = -1.0;
  v[4] = 19.3;
  EXPECT_EQ(quantile(v, 0), -1.);
  EXPECT_EQ(quantile(v, 1), 19.3);
  EXPECT_EQ(quantile(v, 0.2), 1.);
}

TEST(MathFunctions, quantileStdVecDouble) {
  test_quantile_double<std::vector<double> >();
}

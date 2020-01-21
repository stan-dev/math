#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, zeros_array) {
  using stan::math::zeros_array;

  std::vector<double> u0 = zeros_array(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    std::vector<double> u = zeros_array(K);
    std::vector<double> v(K);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_FLOAT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, zeros_array_throw) {
  EXPECT_THROW(stan::math::zeros_array(-1), std::domain_error);
}

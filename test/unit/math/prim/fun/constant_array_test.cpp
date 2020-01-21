#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, constant_array) {
  using stan::math::constant_array;

  std::vector<double> u0 = constant_array(0, 1);
  EXPECT_EQ(0, u0.size());

  double c = 3.14;
  for (int K = 1; K < 5; K++) {
    std::vector<double> u = constant_array(K, c);
    std::vector<double> v(K, c);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_FLOAT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, constant_array_throw) {
  EXPECT_THROW(stan::math::constant_array(-1, 3.14), std::domain_error);
}

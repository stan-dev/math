#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, one_hot_array) {
  for (int K = 1; K < 5; K++) {
    for (int k = 1; k <= K; k++) {
      std::vector<double> v(K);
      v[k - 1] = 1;
      EXPECT_STD_VECTOR_FLOAT_EQ(v, stan::math::one_hot_array(K, k));
    }
  }
}

TEST(MathFunctions, one_hot_array_throw) {
  using stan::math::one_hot_array;
  int K = 5;
  int k = 2;

  EXPECT_THROW(one_hot_array(-1, k), std::domain_error);
  EXPECT_THROW(one_hot_array(0, k), std::domain_error);
  EXPECT_THROW(one_hot_array(K, K + 1), std::domain_error);
  EXPECT_THROW(one_hot_array(K, 0), std::domain_error);
  EXPECT_THROW(one_hot_array(K, -1), std::domain_error);
}

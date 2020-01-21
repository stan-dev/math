#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, one_hot_int_array) {
  using stan::math::one_hot_int_array;

  for (int K = 1; K < 5; K++) {
    for (int k = 1; k <= K; k++) {
      std::vector<int> u = one_hot_int_array(K, k);
      std::vector<int> v(K);
      v[k - 1] = 1;

      EXPECT_EQ(u.size(), v.size());
      for (int i = 0; i < K; i++) {
        EXPECT_EQ(u[i], v[i]);
      }
    }
  }
}

TEST(MathFunctions, one_hot_int_array_throw) {
  using stan::math::one_hot_int_array;
  int K = 5;
  int k = 2;

  EXPECT_THROW(one_hot_int_array(-1, k), std::domain_error);
  EXPECT_THROW(one_hot_int_array(0, k), std::domain_error);
  EXPECT_THROW(one_hot_int_array(K, K + 1), std::domain_error);
  EXPECT_THROW(one_hot_int_array(K, 0), std::domain_error);
  EXPECT_THROW(one_hot_int_array(K, -1), std::domain_error);
}

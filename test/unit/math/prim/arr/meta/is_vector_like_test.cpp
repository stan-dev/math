#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(is_vector_like, vector) {
  using stan::is_vector_like;
  EXPECT_TRUE(is_vector_like<std::vector<double> >::value);
  EXPECT_TRUE(is_vector_like<std::vector<int> >::value);
  EXPECT_TRUE(is_vector_like<const std::vector<double>>::value);
  EXPECT_TRUE(is_vector_like<const std::vector<int> >::value);
  EXPECT_TRUE(is_vector_like<const std::vector<double>&>::value);
  EXPECT_TRUE(is_vector_like<const std::vector<int>&>::value);
}

TEST(is_vector_like, pointer) {
  using stan::is_vector_like;
  EXPECT_TRUE(is_vector_like<double*>::value);
}

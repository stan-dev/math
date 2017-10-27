#include <gtest/gtest.h>
#include <stan/math/prim/arr.hpp>
#include <vector>

TEST(is_vector_like, vector) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<double> >::value);
}

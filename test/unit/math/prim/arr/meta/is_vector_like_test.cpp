#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(is_vector_like, vector) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<double> >::value);
}

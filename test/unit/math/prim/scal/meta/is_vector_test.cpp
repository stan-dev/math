#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_vector) {
  using stan::is_vector;
  EXPECT_FALSE(is_vector<double>::value);
  EXPECT_FALSE(is_vector<int>::value);
  EXPECT_FALSE(is_vector<size_t>::value);

  EXPECT_TRUE(is_vector<std::vector<double>>::value);
  EXPECT_TRUE(is_vector<std::vector<int>>::value);
  EXPECT_TRUE(is_vector<std::vector<double>>::value);
}

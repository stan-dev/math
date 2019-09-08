#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_vector) {
  using stan::is_vector;
  using std::vector;

  EXPECT_TRUE(is_vector<std::vector<double>>::value);
  EXPECT_TRUE(is_vector<std::vector<int>>::value);
  EXPECT_TRUE(is_vector<const std::vector<double>>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>>::value);
  EXPECT_TRUE(is_vector<const std::vector<double>>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>&>::value);
  EXPECT_TRUE(is_vector<std::vector<int>&&>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>&&>::value);
  EXPECT_FALSE(is_vector<double>::value);
}

#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsPrimScal, is_vector) {
  using stan::is_vector;
  EXPECT_FALSE(is_vector<double>::value);
  EXPECT_FALSE(is_vector<int>::value);
  EXPECT_FALSE(is_vector<size_t>::value);

  EXPECT_FALSE(is_vector<const double>::value);
  EXPECT_FALSE(is_vector<const int>::value);
  EXPECT_FALSE(is_vector<const size_t>::value);
  EXPECT_FALSE(is_vector<size_t*>::value);
}

TEST(MetaTraitsPrimArr, is_vector) {
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

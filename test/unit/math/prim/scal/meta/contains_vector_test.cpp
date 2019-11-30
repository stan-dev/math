#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsPrimArr, contains_vector) {
  using stan::contains_vector;
  using std::vector;

  EXPECT_TRUE(contains_vector<std::vector<double> >::value);
  EXPECT_TRUE(contains_vector<std::vector<int> >::value);
  EXPECT_TRUE(contains_vector<const std::vector<double> >::value);
  EXPECT_TRUE(contains_vector<const std::vector<int> >::value);
}

TEST(MetaTraitsPrimScal, contains_vector) {
  using stan::contains_vector;
  EXPECT_FALSE(contains_vector<double>::value);
  EXPECT_FALSE(contains_vector<int>::value);
  EXPECT_FALSE(contains_vector<size_t>::value);

  EXPECT_FALSE(contains_vector<const double>::value);
  EXPECT_FALSE(contains_vector<const int>::value);
  EXPECT_FALSE(contains_vector<const size_t>::value);
}

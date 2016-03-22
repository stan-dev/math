#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, contains_vector) {
  using stan::contains_vector;
  using std::vector;

  EXPECT_TRUE(contains_vector<std::vector<double> >::value);
  EXPECT_TRUE(contains_vector<std::vector<int> >::value);
  EXPECT_TRUE(contains_vector<std::vector<const double> >::value);
  EXPECT_TRUE(contains_vector<std::vector<const int> >::value);
}

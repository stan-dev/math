#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_vector) {
  using stan::is_vector;
  using std::vector;

  EXPECT_TRUE(is_vector<std::vector<double> >::value);
  EXPECT_TRUE(is_vector<std::vector<int> >::value);
  EXPECT_TRUE(is_vector<std::vector<const double> >::value);
  EXPECT_TRUE(is_vector<std::vector<const int> >::value);
  bool temp = is_vector<std::vector<const int>, std::vector<const int>, std::vector<const double> >::value;
  EXPECT_TRUE(temp);
  temp = is_vector<std::vector<const int>, double, double >::value;
  EXPECT_TRUE(temp);
}

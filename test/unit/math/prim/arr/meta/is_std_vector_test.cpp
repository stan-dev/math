#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, contains_vector_true) {
  using stan::is_std_vector;
  using std::vector;

  EXPECT_TRUE(is_std_vector<std::vector<double> >::value);
  EXPECT_TRUE(is_std_vector<std::vector<int> >::value);
  EXPECT_TRUE(is_std_vector<std::vector<const double> >::value);
  EXPECT_TRUE(is_std_vector<std::vector<const int> >::value);
}

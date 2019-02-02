#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, contains_std_vector_false) {
  using stan::contains_std_vector;

  EXPECT_FALSE(contains_std_vector<double>::value);
  EXPECT_FALSE(contains_std_vector<int>::value);
  EXPECT_FALSE(contains_std_vector<const double>::value);
  EXPECT_FALSE(contains_std_vector<const int>::value);
}

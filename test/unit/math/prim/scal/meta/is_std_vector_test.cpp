#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_std_vector_false) {
  using stan::is_std_vector;

  EXPECT_FALSE(is_std_vector<double>::value);
  EXPECT_FALSE(is_std_vector<int>::value);
  EXPECT_FALSE(is_std_vector<const double>::value);
  EXPECT_FALSE(is_std_vector<const int>::value);
}

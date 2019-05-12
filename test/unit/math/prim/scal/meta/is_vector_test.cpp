#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_vector) {
  using stan::is_vector;
  EXPECT_FALSE(is_vector<double>::value);
  EXPECT_FALSE(is_vector<int>::value);
  EXPECT_FALSE(is_vector<size_t>::value);

  EXPECT_FALSE(is_vector<const double>::value);
  EXPECT_FALSE(is_vector<const int>::value);
  EXPECT_FALSE(is_vector<const size_t>::value);

  bool temp = is_vector<double, int, double>::value;
  EXPECT_FALSE(temp);
  temp = is_vector<double, int, const size_t, const size_t, const size_t>::value;
  EXPECT_FALSE(temp);
}

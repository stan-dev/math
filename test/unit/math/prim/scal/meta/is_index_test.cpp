#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_index) {
  using stan::is_index;
  EXPECT_TRUE(is_index<int>::value);
  EXPECT_TRUE(is_index<unsigned int>::value);
  EXPECT_FALSE(is_index<double>::value);
}

#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsPrimScal, is_var) {
  using stan::is_var;
  EXPECT_FALSE(is_var<int>::value);
  EXPECT_FALSE(is_var<double>::value);
}

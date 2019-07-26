#include <stan/math/prim/meta/is_var.hpp>
#include <gtest/gtest.h>


TEST(MetaTraits_scal, is_var) {
  using stan::is_var;
  EXPECT_FALSE(is_var<int>::value);
  EXPECT_FALSE(is_var<double>::value);
}

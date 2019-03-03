


#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<double>::value);
}

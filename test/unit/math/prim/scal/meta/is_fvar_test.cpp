#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<int>::value);
  EXPECT_FALSE(is_fvar<double>::value);
}

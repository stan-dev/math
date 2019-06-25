
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>




TEST(MetaTraits_scal, is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<int>::value);
  EXPECT_FALSE(is_fvar<double>::value);
}

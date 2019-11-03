#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_scalar_container_check) {
  using stan::is_container;
  EXPECT_FALSE(is_container<int>::value);
  EXPECT_FALSE(is_container<double>::value);
}

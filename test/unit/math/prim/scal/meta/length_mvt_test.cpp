#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, length_mvt) {
  using stan::length_mvt;

  double x1;
  EXPECT_THROW(length_mvt(x1), std::invalid_argument);

  int x2;
  EXPECT_THROW(length_mvt(x2), std::invalid_argument);
}

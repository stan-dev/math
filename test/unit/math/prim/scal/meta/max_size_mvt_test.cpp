#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, max_size_mvt) {
  using stan::max_size_mvt;

  double x1 = 0, x2 = 1.0, x3 = -1.0, x4 = 4.0;
  EXPECT_THROW(max_size_mvt(x1, x2), std::invalid_argument);
  EXPECT_THROW(max_size_mvt(x1, x2, x3), std::invalid_argument);
  EXPECT_THROW(max_size_mvt(x1, x2, x3, x4), std::invalid_argument);
}

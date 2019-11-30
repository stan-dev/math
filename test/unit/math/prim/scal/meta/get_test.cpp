#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsPrimArr, get) {
  using stan::get;

  std::vector<double> x(3);
  x[1] = 5.0;
  EXPECT_EQ(5.0, get(x, 1));
}


TEST(MetaTraitsPrimScal, get) {
  using stan::get;

  EXPECT_FLOAT_EQ(2.0, get(2.0, 1));
}

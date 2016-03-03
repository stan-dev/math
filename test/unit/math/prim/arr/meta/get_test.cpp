#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, get) {
  using stan::get;

  std::vector<double> x(3);
  x[1] = 5.0;
  EXPECT_EQ(5.0, get(x,1));
}

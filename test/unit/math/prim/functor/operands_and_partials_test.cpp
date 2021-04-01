#include <stan/math/prim/functor.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, OperandsAndPartials) {
  using stan::math::operands_and_partials;

  operands_and_partials<double> o1(1.0);
  operands_and_partials<double, double, double, double> o2(2.0, 3.0, 4.0, 5.0);

  // This is size 10 because of the two empty broadcast arrays in each edge
  EXPECT_EQ(10, sizeof(o2));

  EXPECT_FLOAT_EQ(27.1, o1.build(27.1));
}

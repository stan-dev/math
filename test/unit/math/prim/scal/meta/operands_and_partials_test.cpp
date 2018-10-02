#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartials) {
  using stan::math::operands_and_partials;

  operands_and_partials<double> o1(1.0);
  operands_and_partials<double, double, double, double> o2(2.0, 3.0, 4.0, 5.0);

  // TODO(Sean): Learn why it's size 6
  EXPECT_GT(sizeof(&o2), sizeof(o2));

  EXPECT_FLOAT_EQ(27.1, o1.build(27.1));
}

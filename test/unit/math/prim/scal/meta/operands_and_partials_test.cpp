#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartials) {
  using stan::math::operands_and_partials;

  operands_and_partials<double> o1(0.0);
  operands_and_partials<double,double,double,double> o2(0.0, 0.0, 0.0, 0.0);

  // TODO(Sean): Learn why it's size 6
  EXPECT_GT(sizeof(&o2), sizeof(o2));

  SUCCEED() << "Construction should be ok.";
}

#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartials) {
  using stan::math::OperandsAndPartials;

  OperandsAndPartials<double> o1;
  OperandsAndPartials<double,double,double,double> o2;

  // TODO(Sean): Learn why it's size 6
  EXPECT_GT(sizeof(&o2), sizeof(o2));

  SUCCEED() << "Construction should be ok.";
}

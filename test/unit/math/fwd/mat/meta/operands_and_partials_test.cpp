#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartialsFvarVec) {
  using stan::math::detail::operands_and_partials;
  using stan::math::fvar;

  std::vector<fvar<double> > x1;
  x1.push_back(fvar<double>(2.0,2.0));
  x1.push_back(fvar<double>(1.0,3.0));

  fvar<double> x2 = 3.0;
  fvar<double> x3 = 5.0;
  x2.d_ = -1.0;
  x3.d_ = 4.0;

  Eigen::VectorXd dx1(2);
  dx1 << 17.0, 13.0;

  operands_and_partials<std::vector<fvar<double> >,fvar<double>,fvar<double> > o(x1, x2, x3);
  o.increment_dx1_vector(0, dx1);
  o.increment_dx2(0, 19.0);
  o.increment_dx2(0, 19.0);
  o.increment_dx3(0, 23.0);
  o.increment_dx3(0, 23.0);
  fvar<double> y = o.build(-1.0);

  EXPECT_FLOAT_EQ(2*17 + 3*13 - 2*19 + 2*4*23,y.d_);
  EXPECT_FLOAT_EQ(-1,y.val_);
}

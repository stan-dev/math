#include <stan/math/prim/fun.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathPrimFun, max_size) {
  using stan::math::max_size;
  double x1 = 1.0;
  double x2 = 0.0;
  double x3 = 2.0;
  double x4 = -3.0;
  double x5 = -11.0;
  EXPECT_EQ(1, max_size(x1, x2));
  EXPECT_EQ(1, max_size(x1, x2, x3));
  EXPECT_EQ(1, max_size(x1, x2, x3, x4));
  EXPECT_EQ(1, max_size(x1, x2, x3, x4, x5));
  std::vector<double> sv{1, 2, 3};
  Eigen::VectorXd ev = Eigen::VectorXd::Zero(4);
  EXPECT_EQ(3, max_size(x1, sv, x3));
  EXPECT_EQ(4, max_size(x1, sv, ev, x3));
  EXPECT_EQ(4, max_size(ev, ev, ev, ev, ev));
  EXPECT_EQ(4, max_size(ev, ev, ev, ev, ev, ev));
  EXPECT_EQ(4, max_size(ev));
  EXPECT_EQ(3, max_size(sv));
  EXPECT_EQ(1, max_size(x1));
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(3, 3);
  EXPECT_EQ(9, max_size(m, x1));
}

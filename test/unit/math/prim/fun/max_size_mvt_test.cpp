#include <stan/math/prim/fun.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathPrimFun, max_size_mvt) {
  using stan::math::max_size_mvt;
  double x1 = 0;
  double x2 = 1.0;
  double x3 = -1.0;
  double x4 = 4.0;
  EXPECT_THROW(max_size_mvt(x1, x2), std::invalid_argument);
  EXPECT_THROW(max_size_mvt(x1, x2, x3), std::invalid_argument);
  EXPECT_THROW(max_size_mvt(x1, x2, x3, x4), std::invalid_argument);
  Eigen::VectorXd ev = Eigen::VectorXd::Zero(3);
  std::vector<Eigen::VectorXd> v0{};
  std::vector<Eigen::VectorXd> v1{ev};
  std::vector<Eigen::VectorXd> v2{ev, ev};
  std::vector<Eigen::VectorXd> v3{ev, ev, ev};
  std::vector<Eigen::VectorXd> v4{ev, ev, ev, ev};
  EXPECT_EQ(0, max_size_mvt(v0));
  EXPECT_EQ(1, max_size_mvt(v1));
  EXPECT_EQ(2, max_size_mvt(v2));
  EXPECT_EQ(3, max_size_mvt(v3));
  EXPECT_EQ(4, max_size_mvt(v4));
  EXPECT_EQ(4, max_size_mvt(v4, v2, v0));
  EXPECT_EQ(2, max_size_mvt(v0, v2, v0));
  EXPECT_EQ(0, max_size_mvt(v0, v0, v0));
  EXPECT_EQ(4, max_size_mvt(v4, v4, v4));
  EXPECT_EQ(4, max_size_mvt(v4, v4));
  EXPECT_EQ(3, max_size_mvt(v2, v3, v2, v3));
  EXPECT_EQ(3, max_size_mvt(v2, v3, v2, v3, v0));
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(3, 3);
  std::vector<Eigen::MatrixXd> vm{m, m};
  EXPECT_EQ(2, max_size_mvt(vm, vm));
  EXPECT_EQ(4, max_size_mvt(v3, m, m, v4));
}

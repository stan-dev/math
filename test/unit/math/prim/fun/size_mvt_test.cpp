#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathPrimFun, size_mvt_scalar) {
  using stan::math::size_mvt;

  double x1;
  EXPECT_THROW(size_mvt(x1), std::invalid_argument);

  int x2;
  EXPECT_THROW(size_mvt(x2), std::invalid_argument);
}

TEST(MathPrimFun, size_mvt_matrices_vectors) {
  using stan::math::size_mvt;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x1(2, 3);
  EXPECT_EQ(1U, size_mvt(x1));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x2(3);
  EXPECT_EQ(3U, size_mvt(x2));

  std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> x3(7);
  EXPECT_EQ(7U, size_mvt(x3));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> x4(7);
  EXPECT_EQ(7U, size_mvt(x4));
}

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, size_mvt_scalar) {
  using stan::math::size_mvt;

  stan::math::var x1;
  EXPECT_THROW(size_mvt(x1), std::invalid_argument);
}

TEST(AgradRev, size_mvt_matrices_vectors) {
  using stan::math::size_mvt;

  stan::math::var_value<Eigen::MatrixXd> x1 = Eigen::MatrixXd(2, 3);
  EXPECT_EQ(1U, size_mvt(x1));

  std::vector<stan::math::var_value<Eigen::MatrixXd>> x2(3);
  EXPECT_EQ(3U, size_mvt(x2));

  std::vector<stan::math::var_value<Eigen::RowVectorXd>> x3(7);
  EXPECT_EQ(7U, size_mvt(x3));

  std::vector<stan::math::var_value<Eigen::VectorXd>> x4(7);
  EXPECT_EQ(7U, size_mvt(x4));
}

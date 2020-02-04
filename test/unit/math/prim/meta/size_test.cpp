#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, size_scalar) {
  using stan::math::size;
  EXPECT_EQ(1U, size(27.0));
  EXPECT_EQ(1U, size(3));
}

TEST(MathMetaPrim, size_vector) {
  using stan::math::size;

  std::vector<int> a(5);
  EXPECT_EQ(5U, size(a));

  std::vector<double> x(10);
  EXPECT_EQ(10U, size(x));

  std::vector<Eigen::MatrixXd> x2(3);
  EXPECT_EQ(3U, size(x2));

  std::vector<Eigen::RowVectorXd> x3(7);
  EXPECT_EQ(7U, size(x3));

  std::vector<Eigen::VectorXd> x4(9);
  EXPECT_EQ(9U, size(x4));
}

TEST(MathMetaPrim, size_matrices) {
  using stan::math::size;

  Eigen::MatrixXd m(2, 3);
  EXPECT_EQ(6U, size(m));

  Eigen::RowVectorXd rv(2);
  EXPECT_EQ(2U, size(rv));

  Eigen::VectorXd v(2);
  EXPECT_EQ(2U, size(v));
}

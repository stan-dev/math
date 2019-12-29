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
  std::vector<double> x(10);
  EXPECT_EQ(10U, size(x));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x2(3);
  EXPECT_EQ(3U, size(x2));

  std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> x3(7);
  EXPECT_EQ(7U, size(x3));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> x4(9);
  EXPECT_EQ(9U, size(x4));
}

TEST(MathMetaPrim, size_matrices) {
  using stan::math::size;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(6U, size(m));

  Eigen::Matrix<double, Eigen::Dynamic, 1> rv(2);
  rv << 1, 2;
  EXPECT_EQ(2U, size(rv));

  Eigen::Matrix<double, 1, Eigen::Dynamic> v(2);
  v << 1, 2;
  EXPECT_EQ(2U, size(v));
}

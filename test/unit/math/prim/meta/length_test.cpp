#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::length;

TEST(MathMetaPrim, length_scalar) {
  using stan::length;
  EXPECT_EQ(1U, size(27.0));
  EXPECT_EQ(1U, size(3));
}

TEST(MathMetaPrim, length_vector) {
  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, size(x));
}

TEST(MathMetaPrim, length_matrices) {
  using stan::length;

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

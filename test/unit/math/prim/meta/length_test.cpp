
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>




TEST(MetaTraits_mat, length) {
  using stan::length;

  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, length(x));

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(6U, length(m));

  Eigen::Matrix<double, Eigen::Dynamic, 1> rv(2);
  rv << 1, 2;
  EXPECT_EQ(2U, length(rv));

  Eigen::Matrix<double, 1, Eigen::Dynamic> v(2);
  v << 1, 2;
  EXPECT_EQ(2U, length(v));

  using stan::length;
  EXPECT_EQ(1U, length(27.0));
  EXPECT_EQ(1U, length(3));
}

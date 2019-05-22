#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>

TEST(MetaTraits, size_of) {
  using stan::size_of;

  std::vector<double> x2(3);
  EXPECT_EQ(3U, size_of(x2));

  std::vector<std::vector<double>> x3(6);
  EXPECT_EQ(6U, size_of(x3));

  std::vector<Eigen::MatrixXd> x4(8);
  EXPECT_EQ(8U, size_of(x4));
}

#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, length_mvt) {
  using stan::length_mvt;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x1(2, 3);
  EXPECT_EQ(1U, length_mvt(x1));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x2(3);
  EXPECT_EQ(3U, length_mvt(x2));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x3(7);
  EXPECT_EQ(7U, length_mvt(x3));
}

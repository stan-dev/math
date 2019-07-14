#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

#include <vector>

TEST(MetaTraits, length_mvt) {
  using stan::length_mvt;

  double x1;
  EXPECT_THROW(length_mvt(x1), std::invalid_argument);

  int x2;
  EXPECT_THROW(length_mvt(x2), std::invalid_argument);
}

TEST(MetaTraits_mat, length_mvt) {
  using stan::length_mvt;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x1(2, 3);
  EXPECT_EQ(1U, length_mvt(x1));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x2(3);
  EXPECT_EQ(3U, length_mvt(x2));

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x3(7);
  EXPECT_EQ(7U, length_mvt(x3));
}

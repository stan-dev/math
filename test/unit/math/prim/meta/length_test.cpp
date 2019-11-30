#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::length;

TEST(MetaTraitsPrimScal, length) {
  using stan::length;
  EXPECT_EQ(1U, length(27.0));
  EXPECT_EQ(1U, length(3));
}

TEST(MetaTraitsPrimArr, length) {
  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, length(x));
}

TEST(MetaTraitsPrimMat, length) {
  using stan::length;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(6U, length(m));

  Eigen::Matrix<double, Eigen::Dynamic, 1> rv(2);
  rv << 1, 2;
  EXPECT_EQ(2U, length(rv));

  Eigen::Matrix<double, 1, Eigen::Dynamic> v(2);
  v << 1, 2;
  EXPECT_EQ(2U, length(v));
}

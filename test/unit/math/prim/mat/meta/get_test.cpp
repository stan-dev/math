#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, get) {
  using stan::get;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 3);
  m << 1, 3, 5, 2, 4, 6;
  EXPECT_EQ(1.0, get(m, 0));
  EXPECT_EQ(3.0, get(m, 2));
  EXPECT_EQ(6.0, get(m, 5));

  Eigen::Matrix<double, Eigen::Dynamic, 1> rv(2);
  rv << 1, 2;
  EXPECT_EQ(1.0, get(rv, 0));
  EXPECT_EQ(2.0, get(rv, 1));

  Eigen::Matrix<double, 1, Eigen::Dynamic> v(2);
  v << 1, 2;
  EXPECT_EQ(1.0, get(v, 0));
  EXPECT_EQ(2.0, get(v, 1));
}

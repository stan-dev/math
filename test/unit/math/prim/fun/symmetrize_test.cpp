#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(Symmetrize, Error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::symmetrize;

  Matrix<double, Dynamic, Dynamic> m(3, 4);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      m(i, j) = (i + 1) * (j + 1);
  EXPECT_THROW(symmetrize(m), std::invalid_argument);
}

TEST(Symmetrize, Value) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::symmetrize;
  using stan::math::is_symmetric;

  Matrix<double, Dynamic, Dynamic> m(4, 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      m(i, j) = (i + 1) * (j + 1);
  EXPECT_TRUE(is_symmetric(symmetrize(m)));
}


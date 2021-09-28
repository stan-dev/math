#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradPrimMatrix, initialize_fill) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::initialize_fill;
  using std::vector;

  double x;
  double y = 10;
  initialize_fill(x, y);
  EXPECT_FLOAT_EQ(10.0, x);

  std::vector<double> z(2);
  double a = 15;
  initialize_fill(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0]);
  EXPECT_FLOAT_EQ(15.0, z[1]);
  EXPECT_EQ(2U, z.size());

  Matrix<double, Dynamic, Dynamic> m(2, 3);
  initialize_fill(m, static_cast<double>(12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j));

  Matrix<double, Dynamic, 1> rv(3);
  initialize_fill(rv, static_cast<double>(13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i));

  Matrix<double, 1, Dynamic> v(4);
  initialize_fill(v, static_cast<double>(22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i));

  vector<vector<double> > d(3, vector<double>(2));
  initialize_fill(d, static_cast<double>(54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j]);
}
TEST(AgradPrimMatrix, initialize_fillDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::initialize_fill;
  Matrix<double, Dynamic, 1> y(3);
  initialize_fill(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0]);
}

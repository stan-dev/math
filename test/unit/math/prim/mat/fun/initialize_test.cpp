#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrix, initialize) {
  // 2nd template
  using stan::math::initialize;
  double x;
  double y = 10;
  // template 2
  initialize(x, y);
  EXPECT_FLOAT_EQ(y, x);

  int z = 5;
  // template 2
  initialize(y, z);
  EXPECT_FLOAT_EQ(z, y);
}

TEST(MathMatrix, initMatrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::initialize;
  Matrix<double, Dynamic, Dynamic> m(3, 2);
  // template 3, 2
  initialize(m, 13.2);
  for (int i = 0; i < m.size(); ++i)
    EXPECT_FLOAT_EQ(m(i), 13.2);

  Matrix<double, Dynamic, 1> v(3);
  // template 3, 2
  initialize(v, 2);
  for (int i = 0; i < v.size(); ++i)
    EXPECT_FLOAT_EQ(v(i), 2);

  Matrix<double, 1, Dynamic> rv(3);
  // template 3, 2
  initialize(rv, 12);
  for (int i = 0; i < v.size(); ++i)
    EXPECT_FLOAT_EQ(rv(i), 12);
}

TEST(MathMatrix, initStdVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::initialize;
  using std::vector;

  vector<double> x(3);
  initialize(x, 2.2);
  for (size_t i = 0; i < 3; ++i)
    // template 4, 2
    EXPECT_FLOAT_EQ(2.2, x[i]);

  vector<Matrix<double, Dynamic, Dynamic> > z(
      4, Matrix<double, Dynamic, Dynamic>(3, 2));
  initialize(z, 3.7);
  for (size_t i = 0; i < 4; ++i)
    for (int m = 0; m < 3; ++m)
      for (int n = 0; n < 2; ++n)
        EXPECT_FLOAT_EQ(3.7, z[i](m, n));
}

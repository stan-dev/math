#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, value_of_rec) {
  using stan::math::value_of_rec;
  double x = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of_rec(x));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(5));
}

TEST(MathMatrixPrimArr, value_of_rec) {
  using stan::math::value_of_rec;
  using std::vector;

  vector<double> a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(i + 1);

  vector<double> b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(i + 1);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i], d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i], d_a[i]);
}

TEST(MathMatrixPrimMat, value_of_rec) {
  using stan::math::value_of_rec;

  Eigen::Matrix<double, 2, 5> a;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++)
      a(i, j) = i * 5 + j;
  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 1; j++)
      b(i, j) = 10 + i * 5 + j;

  Eigen::MatrixXd d_a = value_of_rec(a);
  Eigen::VectorXd d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b(i), d_b(i));

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j)
      EXPECT_FLOAT_EQ(a(i, j), d_a(i, j));
}

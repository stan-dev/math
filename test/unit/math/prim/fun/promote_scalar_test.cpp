#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MathFunctionsPromoteScalar, Match) {
  using stan::math::promote_scalar;
  EXPECT_FLOAT_EQ(1.3, promote_scalar<double>(1.3));
  EXPECT_EQ(3, promote_scalar<int>(3));

  EXPECT_TYPE(double, promote_scalar<double>(2.3));
  EXPECT_TYPE(int, promote_scalar<int>(2));
}
TEST(MathFunctionsPromoteScalar, Mismatch) {
  using stan::math::promote_scalar;
  EXPECT_FLOAT_EQ(2.0, promote_scalar<double>(2));
  EXPECT_TYPE(double, promote_scalar<double>(2));
}

TEST(MathFunctionsPromoteScalar, VectorMismatch) {
  using stan::math::promote_scalar;
  std::vector<int> x;
  x.push_back(1);
  x.push_back(2);
  std::vector<double> y = promote_scalar<double>(x);
  EXPECT_EQ(2U, y.size());
  EXPECT_FLOAT_EQ(1.0, y[0]);
  EXPECT_FLOAT_EQ(2.0, y[1]);
}
TEST(MathFunctionsPromoteScalar, VectorMatch) {
  using stan::math::promote_scalar;
  std::vector<int> x;
  x.push_back(13);
  x.push_back(27);
  std::vector<int> y = promote_scalar<int>(x);
  EXPECT_EQ(2U, y.size());
  EXPECT_EQ(13, y[0]);
  EXPECT_EQ(27, y[1]);
}
TEST(MathFunctionsPromoteScalar, Vector2Mismatch) {
  using stan::math::promote_scalar;
  using std::vector;
  vector<vector<int> > x(2);
  x[0].push_back(1);
  x[0].push_back(2);
  x[0].push_back(3);
  x[1].push_back(4);
  x[1].push_back(5);
  x[1].push_back(6);

  vector<vector<double> > y = promote_scalar<double>(x);
  EXPECT_EQ(2U, y.size());
  EXPECT_EQ(3, y[0].size());
  EXPECT_FLOAT_EQ(1.0, y[0][0]);
  EXPECT_FLOAT_EQ(2.0, y[0][1]);
  EXPECT_FLOAT_EQ(3.0, y[0][2]);
  EXPECT_FLOAT_EQ(4.0, y[1][0]);
  EXPECT_FLOAT_EQ(5.0, y[1][1]);
  EXPECT_FLOAT_EQ(6.0, y[1][2]);
}
TEST(MathFunctionsPromoteScalar, Vector2Match) {
  using stan::math::promote_scalar;
  using std::vector;
  vector<vector<double> > x(2);
  x[0].push_back(1.1);
  x[0].push_back(2.2);
  x[0].push_back(3.3);
  x[1].push_back(4.4);
  x[1].push_back(5.5);
  x[1].push_back(6.6);

  vector<vector<double> > y = promote_scalar<double>(x);
  EXPECT_EQ(2U, y.size());
  EXPECT_EQ(3, y[0].size());
  EXPECT_FLOAT_EQ(1.1, y[0][0]);
  EXPECT_FLOAT_EQ(2.2, y[0][1]);
  EXPECT_FLOAT_EQ(3.3, y[0][2]);
  EXPECT_FLOAT_EQ(4.4, y[1][0]);
  EXPECT_FLOAT_EQ(5.5, y[1][1]);
  EXPECT_FLOAT_EQ(6.6, y[1][2]);
}

TEST(MathMatrixPromoteScalar, MatrixMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<int, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Matrix<double, Dynamic, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(6, y.size());
  EXPECT_FLOAT_EQ(1.0, y(0, 0));
  EXPECT_FLOAT_EQ(2.0, y(0, 1));
  EXPECT_FLOAT_EQ(3.0, y(0, 2));
  EXPECT_FLOAT_EQ(4.0, y(1, 0));
  EXPECT_FLOAT_EQ(5.0, y(1, 1));
  EXPECT_FLOAT_EQ(6.0, y(1, 2));
}
TEST(MathMatrixPromoteScalar, MatrixMatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<double, Dynamic, Dynamic> x(2, 3);
  x << 1.1, 2.2, 3.3, 4.4, 5.5, 6.6;

  Matrix<double, Dynamic, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(6, y.size());
  EXPECT_FLOAT_EQ(1.1, y(0, 0));
  EXPECT_FLOAT_EQ(2.2, y(0, 1));
  EXPECT_FLOAT_EQ(3.3, y(0, 2));
  EXPECT_FLOAT_EQ(4.4, y(1, 0));
  EXPECT_FLOAT_EQ(5.5, y(1, 1));
  EXPECT_FLOAT_EQ(6.6, y(1, 2));
}
TEST(MathMatrixPromoteScalar, VectorMatrixMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;
  using std::vector;

  vector<Matrix<int, Dynamic, Dynamic> > x(2);
  Matrix<int, Dynamic, Dynamic> x0(2, 3);
  x0 << 1, 2, 3, 4, 5, 6;
  x[0] = x0;
  Matrix<int, Dynamic, Dynamic> x1(2, 3);
  x1 << 10, 20, 30, 40, 50, 60;
  x[1] = x1;

  vector<Matrix<double, Dynamic, Dynamic> > y = promote_scalar<double>(x);
  EXPECT_EQ(2, y.size());

  EXPECT_FLOAT_EQ(1.0, y[0](0, 0));
  EXPECT_FLOAT_EQ(2.0, y[0](0, 1));
  EXPECT_FLOAT_EQ(3.0, y[0](0, 2));
  EXPECT_FLOAT_EQ(4.0, y[0](1, 0));
  EXPECT_FLOAT_EQ(5.0, y[0](1, 1));
  EXPECT_FLOAT_EQ(6.0, y[0](1, 2));
  EXPECT_FLOAT_EQ(10.0, y[1](0, 0));
  EXPECT_FLOAT_EQ(20.0, y[1](0, 1));
  EXPECT_FLOAT_EQ(30.0, y[1](0, 2));
  EXPECT_FLOAT_EQ(40.0, y[1](1, 0));
  EXPECT_FLOAT_EQ(50.0, y[1](1, 1));
  EXPECT_FLOAT_EQ(60.0, y[1](1, 2));
}

TEST(MathMatrixPromoteScalar, ColVectorMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<int, Dynamic, 1> x(3);
  x << 1, 2, 3;

  Matrix<double, Dynamic, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(1.0, y(0));
  EXPECT_FLOAT_EQ(2.0, y(1));
  EXPECT_FLOAT_EQ(3.0, y(2));
}
TEST(MathMatrixPromoteScalar, ColVectorMatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<double, Dynamic, 1> x(3);
  x << 1.1, 2.2, 3.3;

  Matrix<double, Dynamic, 1> y = promote_scalar<double>(x);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(1.1, y(0));
  EXPECT_FLOAT_EQ(2.2, y(1));
  EXPECT_FLOAT_EQ(3.3, y(2));
}
TEST(MathMatrixPromoteScalar, VectorColVectorMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;
  using std::vector;

  vector<Matrix<int, Dynamic, 1> > x(2);
  Matrix<int, Dynamic, 1> x0(3);
  x0 << 1, 2, 3;
  x[0] = x0;
  Matrix<int, Dynamic, 1> x1(3);
  x1 << 10, 20, 30;
  x[1] = x1;

  vector<Matrix<double, Dynamic, 1> > y = promote_scalar<double>(x);
  EXPECT_EQ(2, y.size());

  EXPECT_FLOAT_EQ(1.0, y[0](0));
  EXPECT_FLOAT_EQ(2.0, y[0](1));
  EXPECT_FLOAT_EQ(3.0, y[0](2));
  EXPECT_FLOAT_EQ(10, y[1](0));
  EXPECT_FLOAT_EQ(20, y[1](1));
  EXPECT_FLOAT_EQ(30, y[1](2));
}

TEST(MathMatrixPromoteScalar, RowVectorMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<int, 1, Dynamic> x(3);
  x << 1, 2, 3;

  Matrix<double, Dynamic, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(1.0, y(0));
  EXPECT_FLOAT_EQ(2.0, y(1));
  EXPECT_FLOAT_EQ(3.0, y(2));
}
TEST(MathMatrixPromoteScalar, RowVectorMatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;

  Matrix<double, 1, Dynamic> x(3);
  x << 1.1, 2.2, 3.3;

  Matrix<double, 1, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(1.1, y(0));
  EXPECT_FLOAT_EQ(2.2, y(1));
  EXPECT_FLOAT_EQ(3.3, y(2));
}
TEST(MathMatrixPromoteScalar, VectorRowVectorMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::promote_scalar;
  using std::vector;

  vector<Matrix<int, 1, Dynamic> > x(2);
  Matrix<int, 1, Dynamic> x0(3);
  x0 << 1, 2, 3;
  x[0] = x0;
  Matrix<int, 1, Dynamic> x1(3);
  x1 << 10, 20, 30;
  x[1] = x1;

  vector<Matrix<double, 1, Dynamic> > y = promote_scalar<double>(x);
  EXPECT_EQ(2, y.size());

  EXPECT_FLOAT_EQ(1.0, y[0](0));
  EXPECT_FLOAT_EQ(2.0, y[0](1));
  EXPECT_FLOAT_EQ(3.0, y[0](2));
  EXPECT_FLOAT_EQ(10, y[1](0));
  EXPECT_FLOAT_EQ(20, y[1](1));
  EXPECT_FLOAT_EQ(30, y[1](2));
}

TEST(MathMatrixPromoteScalar, Array) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using stan::math::promote_scalar;

  Array<int, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Array<double, Dynamic, Dynamic> y = promote_scalar<double>(x);
  EXPECT_EQ(6, y.size());
  EXPECT_FLOAT_EQ(1.0, y(0, 0));
  EXPECT_FLOAT_EQ(2.0, y(0, 1));
  EXPECT_FLOAT_EQ(3.0, y(0, 2));
  EXPECT_FLOAT_EQ(4.0, y(1, 0));
  EXPECT_FLOAT_EQ(5.0, y(1, 1));
  EXPECT_FLOAT_EQ(6.0, y(1, 2));
}

TEST(MathMatrixPromoteScalar, Expression) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using stan::math::promote_scalar;

  Array<int, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Array<double, Dynamic, Dynamic> y = promote_scalar<double>(x - 1);
  EXPECT_EQ(6, y.size());
  EXPECT_FLOAT_EQ(0.0, y(0, 0));
  EXPECT_FLOAT_EQ(1.0, y(0, 1));
  EXPECT_FLOAT_EQ(2.0, y(0, 2));
  EXPECT_FLOAT_EQ(3.0, y(1, 0));
  EXPECT_FLOAT_EQ(4.0, y(1, 1));
  EXPECT_FLOAT_EQ(5.0, y(1, 2));
}

#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

template <int N>
void test_print_mat_size(const std::string& expected) {
  using stan::math::print_mat_size;
  std::stringstream ss;
  stan::math::print_mat_size<N>(ss);
  EXPECT_EQ(expected, ss.str());
}

TEST(MathMatrixAssign, print_mat_size) {
  test_print_mat_size<-1>("dynamically sized");
  test_print_mat_size<0>("0");
  test_print_mat_size<1>("1");
  test_print_mat_size<10>("10");
}

TEST(MathMatrixAssign, realToDouble) {
  double a = 2.7;
  int b = 0;
  EXPECT_NO_THROW(stan::math::assign(b, a));
  EXPECT_NO_THROW(stan::math::assign(a, b));
}

TEST(MathMatrixAssign, matSizeMismatch) {
  using Eigen::Matrix;
  using stan::math::assign;
  Matrix<double, 1, -1> x(3);
  Matrix<double, 1, -1> y(3);
  EXPECT_NO_THROW(assign(x, y));
  Matrix<double, 1, -1> z(10);
  EXPECT_THROW_MSG(assign(x, z), std::exception,
                   "Columns of left-hand-side (3) and columns"
                   " of right-hand-side (10) must match in size");

  Matrix<double, -1, 1> a(5);
  Matrix<double, -1, 1> b(6);
  EXPECT_THROW_MSG(assign(a, b), std::exception,
                   "Rows of left-hand-side (5) and rows"
                   " of right-hand-side (6) must match in size");

  std::vector<double> u(5);
  std::vector<double> v(7);
  EXPECT_THROW_MSG(assign(u, v), std::exception,
                   "size of left-hand side (5) and size of"
                   " right-hand side (7) must match in size");

  Matrix<double, -1, -1> m1(4, 5);
  Matrix<double, -1, -1> m2(5, 4);
  EXPECT_THROW_MSG(assign(m1, m2), std::exception,
                   "Rows of left-hand-side (4) and rows of"
                   " right-hand-side (5) must match in size");

  Matrix<double, -1, -1> m3(10, 10);
  Matrix<double, -1, -1> m4(2, 3);
  EXPECT_THROW_MSG(assign(m3.block(1, 1, 7, 3), m4), std::exception,
                   "left-hand side rows (7) and right-hand"
                   " side rows (2) must match in size");

  EXPECT_THROW_MSG(assign(m3.block(1, 1, 2, 5), m4), std::exception,
                   "assign: left-hand side cols (5) and"
                   " right-hand side cols (3) must match in size");
}

TEST(MathMatrixAssign, test_int) {
  using stan::math::assign;
  int a;
  int b = 5;
  assign(a, b);
  EXPECT_EQ(5, a);
  EXPECT_EQ(5, b);

  assign(a, 12);
  EXPECT_EQ(a, 12);
}
TEST(MathMatrixAssign, test_double) {
  using stan::math::assign;

  double a;
  int b = 5;
  double c = 5.0;
  assign(a, b);
  EXPECT_FLOAT_EQ(5.0, a);
  EXPECT_FLOAT_EQ(5.0, b);

  assign(a, c);
  EXPECT_FLOAT_EQ(5.0, a);
  EXPECT_FLOAT_EQ(5.0, b);

  assign(a, 5.2);
  EXPECT_FLOAT_EQ(5.2, a);
}
TEST(MathMatrixAssign, vectorDouble) {
  using stan::math::assign;
  using std::vector;

  vector<double> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;

  vector<double> x(3);
  assign(x, y);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(y[i], x[i]);

  vector<double> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  vector<int> ns(3);
  ns[0] = 1;
  ns[1] = -10;
  ns[2] = 500;

  assign(x, ns);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, ns.size());
  for (size_t i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(ns[i], x[i]);
}

TEST(MathMatrixAssign, eigenRowVectorDoubleToDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<double, 1, Dynamic> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;

  Matrix<double, 1, Dynamic> x(3);
  assign(x, y);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(y[i], x[i]);
}
TEST(MathMatrixAssign, eigenRowVectorIntToDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<double, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<int, 1, Dynamic> ns(3);
  ns[0] = 1;
  ns[1] = -10;
  ns[2] = 500;

  assign(x, ns);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, ns.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(ns[i], x[i]);
}
TEST(MathMatrixAssign, eigenRowVectorShapeMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<double, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<double, 1, Dynamic> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  Matrix<double, Dynamic, 1> zz(3);
  zz << 1, 2, 3;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);

  Matrix<double, Dynamic, Dynamic> zzz(3, 1);
  zzz << 1, 2, 3;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);

  Matrix<double, Dynamic, Dynamic> zzzz(1, 3);
  EXPECT_THROW(assign(x, zzzz), std::invalid_argument);
}

TEST(MathMatrixAssign, eigenMatrixDoubleToDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<double, Dynamic, Dynamic> y(3, 2);
  y << 1.2, 100, -5.1, 12, 1000, -5100;

  Matrix<double, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(y(i), x(i));
}
TEST(MathMatrixAssign, eigenMatrixIntToDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<int, Dynamic, Dynamic> y(3, 2);
  y << 1, 2, 3, 4, 5, 6;

  Matrix<double, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(y(i), x(i));
}
TEST(MathMatrixAssign, eigenMatrixShapeMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;

  Matrix<double, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Matrix<double, 1, Dynamic> z(6);
  z << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, z), std::invalid_argument);
  EXPECT_THROW(assign(z, x), std::invalid_argument);

  Matrix<double, Dynamic, 1> zz(6);
  zz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);
  EXPECT_THROW(assign(zz, x), std::invalid_argument);

  Matrix<double, Dynamic, Dynamic> zzz(6, 1);
  zzz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);
  EXPECT_THROW(assign(zzz, x), std::invalid_argument);
}

TEST(MathMatrixPrimMat, block) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;
  using stan::math::get_base1_lhs;

  Matrix<double, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<double, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0));
  EXPECT_FLOAT_EQ(100.0, m(0, 1));
  EXPECT_FLOAT_EQ(1000.0, m(0, 2));
}

TEST(MathMatrixPrimMat, block2) {
  using Eigen::MatrixXd;
  using stan::math::assign;

  MatrixXd a(2, 3);
  a << 1, 2, 3, 4, 5, 6;

  MatrixXd b(2, 2);
  b << 10, 20, 30, 40;

  assign(a.block(0, 0, 2, 2), b);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(a(i, j), b(i, j));

  EXPECT_FLOAT_EQ(a(0, 2), 3.0);
  EXPECT_FLOAT_EQ(a(1, 2), 6.0);
}

TEST(MathMatrixPrimMat, vectorVector) {
  using stan::math::assign;
  using std::vector;
  vector<vector<double> > x(3, vector<double>(2));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      x[i][j] = (i + 1) * (j - 10);

  vector<vector<double> > y(3, vector<double>(2));

  assign(y, x);
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3U; ++i) {
    EXPECT_EQ(2U, y[i].size());
    for (size_t j = 0; j < 2U; ++j) {
      EXPECT_FLOAT_EQ(x[i][j], y[i][j]);
    }
  }
}

TEST(MathMatrixPrimMat, vectorVectorVector) {
  using stan::math::assign;
  using std::vector;
  vector<vector<vector<double> > > x(
      4, vector<vector<double> >(3, vector<double>(2)));
  for (size_t k = 0; k < 4; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j)
        x[k][i][j] = (i + 1) * (j - 10) * (20 * k + 100);

  vector<vector<vector<double> > > y(
      4, vector<vector<double> >(3, vector<double>(2)));

  assign(y, x);
  EXPECT_EQ(4U, y.size());
  for (size_t k = 0; k < 4U; ++k) {
    EXPECT_EQ(3U, y[k].size());
    for (size_t i = 0; i < 3U; ++i) {
      EXPECT_EQ(2U, y[k][i].size());
      for (size_t j = 0; j < 2U; ++j) {
        EXPECT_FLOAT_EQ(x[k][i][j], y[k][i][j]);
      }
    }
  }
}

TEST(MathMatrixPrimMat, vectorEigenVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;
  using std::vector;

  vector<Matrix<double, Dynamic, 1> > x(2, Matrix<double, Dynamic, 1>(3));
  for (size_t i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      x[i](j) = (i + 1) * (10 * j + 2);
  vector<Matrix<double, Dynamic, 1> > y(2, Matrix<double, Dynamic, 1>(3));

  assign(y, x);

  EXPECT_EQ(2U, y.size());
  for (size_t i = 0; i < 2U; ++i) {
    EXPECT_EQ(3U, y[i].size());
    for (size_t j = 0; j < 3U; ++j) {
      EXPECT_FLOAT_EQ(x[i](j), y[i](j));
    }
  }
}

TEST(MathMatrixPrimMat, getAssignRow) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::assign;
  using stan::math::get_base1_lhs;

  Matrix<double, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<double, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0));
  EXPECT_FLOAT_EQ(100.0, m(0, 1));
  EXPECT_FLOAT_EQ(1000.0, m(0, 2));
}

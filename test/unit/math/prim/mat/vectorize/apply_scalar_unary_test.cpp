#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

TEST(MathPrimMatVectorize, applyScalarUnaryReturnT) {
  using stan::test::expect_match_return_t;
  expect_match_return_t<double, int>();
  typedef stan::math::apply_scalar_unary<stan::math::foo_fun,
                                         std::vector<int> >::return_t
    vec_double_return_t;
  vec_double_return_t f;
  f.push_back(3.7);
  EXPECT_FLOAT_EQ(3.7, f[0]);
  expect_match_return_t<std::vector<double>, std::vector<int> >();
}

TEST(MathPrimMatVectorize, applyScalarUnaryScalar) {
  using stan::math::foo;

  double exp_3 = foo(3);
  EXPECT_FLOAT_EQ(std::exp(3.0), exp_3);

  double exp_4_d = foo(4.0);
  EXPECT_FLOAT_EQ(std::exp(4.0), exp_4_d);

  EXPECT_FLOAT_EQ(std::exp(2.3), foo(2.3));
}
TEST(MathPrimMatVectorize, applyScalarUnaryStdVector) {
  using stan::math::foo;
  using std::vector;

  vector<double> y;
  y.push_back(1.3);
  y.push_back(-2.9);

  vector<double> fy = foo(y);
  EXPECT_EQ(2, fy.size());
  EXPECT_FLOAT_EQ(std::exp(y[0]), fy[0]);
  EXPECT_FLOAT_EQ(std::exp(y[1]), fy[1]);

  vector<double> z1;
  z1.push_back(1.1);
  z1.push_back(1.2);
  vector<double> z2;
  z2.push_back(2.1);
  z2.push_back(2.2);
  vector<vector<double> > z;
  z.push_back(z1);
  z.push_back(z2);
  vector<vector<double> > fz = foo(z);
  EXPECT_EQ(z.size(), fz.size());
  for (size_t i = 0; i < fz.size(); ++i) {
    EXPECT_EQ(z[i].size(), fz[i].size());
    for (size_t j = 0; j < z[i].size(); ++j)
      EXPECT_FLOAT_EQ(std::exp(z[i][j]), fz[i][j]);
  }

  vector<int> u1;
  u1.push_back(2);
  u1.push_back(3);
  vector<double> u2 = foo(u1);
  EXPECT_EQ(2, u2.size());
  for (int i = 0; i < u2.size(); ++i) {
    EXPECT_FLOAT_EQ(std::exp(u1[i]), u2[i]);
  }
}    
TEST(MathPrimMatVectorize, applyScalarUnaryMatrix) {
  using stan::math::foo;
  using Eigen::MatrixXd;

  MatrixXd a(3,2);
  a << 1, 2, 3, 4, 5, 6;
  MatrixXd fa = foo(a);
  EXPECT_EQ(a.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(a(i)), fa(i));

  MatrixXd fab = foo(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(std::exp(a(1,1)), fab(0,0));
}
TEST(MathPrimMatVectorize, applyScalarUnaryVector) {
  using stan::math::foo;
  using Eigen::VectorXd;
  VectorXd b(6);
  b << 1, 2, 3, 4, 5, 6;
  VectorXd fb = foo(b);
  EXPECT_EQ(b.size(), fb.size());
  for (int i = 0; i < fb.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(b(i)), fb(i));
}
TEST(MathPrimMatVectorize, applyScalarUnaryRowVector) {
  using stan::math::foo;
  using Eigen::RowVectorXd;
  RowVectorXd c(6);
  c << 1, 2, 3, 4, 5, 6;
  RowVectorXd fc = foo(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(c(i)), fc(i));
}
TEST(MathPrimMatVectorize, applyScalarUnaryMixedMatrixStdVector) {
  using stan::math::foo;
  using std::vector;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;

  RowVectorXd d1(3);
  d1 << 1, 2, 3;
  RowVectorXd d2(3);
  d2 << 10, 20, 30;
  vector<RowVectorXd> d;
  d.push_back(d1);
  d.push_back(d2);
  vector<RowVectorXd> fd = foo(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(std::exp(d[i](j)), fd[i](j));
  }
}



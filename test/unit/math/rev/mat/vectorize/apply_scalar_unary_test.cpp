#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

TEST(RevMatVectorize,applyScalarUnary) {
  using stan::math::var;
  using stan::math::foo;
  using stan::test::expect_match_return_t;
  var three_var = 3;
  var exp_3_v = foo(three_var);
  EXPECT_FLOAT_EQ(std::exp(3.0), exp_3_v.val());

  expect_match_return_t<var, var>();
  expect_match_return_t<std::vector<var>, std::vector<var> >();
}

TEST(MathPrimMatVectorize, applyScalarUnaryStdVector) {
  using stan::math::foo;
  using std::vector;
  using stan::math::var;
  using stan::math::foo;

  vector<var> y;
  y.push_back(1.3);
  y.push_back(-2.9);

  vector<var> fy = foo(y);
  EXPECT_EQ(2, fy.size());
  EXPECT_FLOAT_EQ(std::exp(y[0].val()), fy[0].val());
  EXPECT_FLOAT_EQ(std::exp(y[1].val()), fy[1].val());

  vector<var> z1;
  z1.push_back(1.1);
  z1.push_back(1.2);
  vector<var> z2;
  z2.push_back(2.1);
  z2.push_back(2.2);
  vector<vector<var> > z;
  z.push_back(z1);
  z.push_back(z2);
  vector<vector<var> > fz = foo(z);
  EXPECT_EQ(z.size(), fz.size());
  for (size_t i = 0; i < fz.size(); ++i) {
    EXPECT_EQ(z[i].size(), fz[i].size());
    for (size_t j = 0; j < z[i].size(); ++j)
      EXPECT_FLOAT_EQ(std::exp(z[i][j].val()), fz[i][j].val());
  }
}    
TEST(MathPrimMatVectorize, applyScalarUnaryMatrix) {
  using stan::math::var;
  using stan::math::foo;
  typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;

  MatrixXvar a(3,2);
  a << 1, 2, 3, 4, 5, 6;
  MatrixXvar fa = foo(a);
  EXPECT_EQ(a.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(a(i).val()), fa(i).val());

  MatrixXvar fab = foo(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(std::exp(a(1,1).val()), fab(0,0).val());
}
TEST(MathPrimMatVectorize, applyScalarUnaryVector) {
  using stan::math::var;
  using stan::math::foo;
  typedef Eigen::Matrix<var, Eigen::Dynamic, 1> VectorXvar;

  VectorXvar b(6);
  b << 1, 2, 3, 4, 5, 6;
  VectorXvar fb = foo(b);
  EXPECT_EQ(b.size(), fb.size());
  for (int i = 0; i < fb.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(b(i).val()), fb(i).val());
}
TEST(MathPrimMatVectorize, applyScalarUnaryRowVector) {
  using stan::math::var;
  using stan::math::foo;
  typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RowVectorXvar;

  RowVectorXvar c(6);
  c << 1, 2, 3, 4, 5, 6;
  RowVectorXvar fc = foo(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i)
    EXPECT_FLOAT_EQ(std::exp(c(i).val()), fc(i).val());
}
TEST(MathPrimMatVectorize, applyScalarUnaryMixedMatrixStdVector) {
  using stan::math::var;
  using stan::math::foo;
  using std::vector;
  typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RowVectorXvar;

  RowVectorXvar d1(3);
  d1 << 1, 2, 3;
  RowVectorXvar d2(3);
  d2 << 10, 20, 30;
  vector<RowVectorXvar> d;
  d.push_back(d1);
  d.push_back(d2);
  vector<RowVectorXvar> fd = foo(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(std::exp(d[i](j).val()), fd[i](j).val());
  }
}

#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

// TODO(carpenter): figure out some way to deal with constraints (transforms)
//                  including testing out of support values
// ?? wrap this all up in a general transform framework



template <typename F>
void expect_scalar_value() {
  double v = F::template apply<double>(3);
  EXPECT_FLOAT_EQ(F::apply_base(3.0), v);

  double v2 = F::template apply<double>(4.0);
  EXPECT_FLOAT_EQ(F::apply_base(4.0), v2);

  EXPECT_FLOAT_EQ(F::apply_base(2.3), F::template apply<double>(2.3));
}

template <typename F>
void expect_std_vectors_value() {
    using stan::math::foo;
  using std::vector;

  vector<double> y;
  y.push_back(1.3);
  y.push_back(-2.9);

  vector<double> fy = F::template apply<vector<double> >(y);
  EXPECT_EQ(2, fy.size());
  EXPECT_FLOAT_EQ(F::apply_base(y[0]), fy[0]);
  EXPECT_FLOAT_EQ(F::apply_base(y[1]), fy[1]);

  vector<double> z1;
  z1.push_back(1.1);
  z1.push_back(1.2);
  vector<double> z2;
  z2.push_back(2.1);
  z2.push_back(2.2);
  vector<vector<double> > z;
  z.push_back(z1);
  z.push_back(z2);
  vector<vector<double> > fz = F::template apply<vector<vector<double> > >(z);
  EXPECT_EQ(z.size(), fz.size());
  for (size_t i = 0; i < fz.size(); ++i) {
    EXPECT_EQ(z[i].size(), fz[i].size());
    for (size_t j = 0; j < z[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(z[i][j]), fz[i][j]);
  }

  vector<int> u1;
  u1.push_back(2);
  u1.push_back(3);
  vector<double> u2 = F::template apply<vector<double> >(u1);
  EXPECT_EQ(2, u2.size());
  for (int i = 0; i < u2.size(); ++i) {
    EXPECT_FLOAT_EQ(F::apply_base(u1[i]), u2[i]);
  }
}

template <typename F>
void expect_matrix_value() {
  using Eigen::MatrixXd;

  MatrixXd a(3,2);
  a << 1, 2, 3, 4, 5, 6;
  MatrixXd fa = F::template apply<MatrixXd>(a);
  EXPECT_EQ(a.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(a(i)), fa(i));
  
  MatrixXd fab = F::template apply<MatrixXd>(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(std::exp(a(1,1)), fab(0,0));
}

template <typename F>
void expect_vector_value() {
  using Eigen::VectorXd;

  VectorXd b(6);
  b << 1, 2, 3, 4, 5, 6;
  VectorXd fb = F::template apply<VectorXd>(b);
  EXPECT_EQ(b.size(), fb.size());
  for (int i = 0; i < fb.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(b(i)), fb(i));
}

template <typename F>
void expect_row_vector_value() {
  using Eigen::RowVectorXd;

  RowVectorXd c(6);
  c << 1, 2, 3, 4, 5, 6;
  RowVectorXd fc = F::template apply<RowVectorXd>(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(c(i)), fc(i));
}

template <typename F>
void expect_std_vector_row_vector_value() {
  using std::vector;
  using Eigen::RowVectorXd;

  RowVectorXd d1(3);
  d1 << 1, 2, 3;
  RowVectorXd d2(3);
  d2 << 10, 20, 30;
  vector<RowVectorXd> d;
  d.push_back(d1);
  d.push_back(d2);
  vector<RowVectorXd> fd = F::template apply<vector<RowVectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
  }
}

template <typename F>
void expect_std_vector_vector_value() {
  using std::vector;
  using Eigen::VectorXd;

  VectorXd d1(3);
  d1 << 1, 2, 3;
  VectorXd d2(3);
  d2 << 10, 20, 30;
  vector<VectorXd> d;
  d.push_back(d1);
  d.push_back(d2);
  vector<VectorXd> fd = F::template apply<vector<VectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
  }
}

template <typename F>
void expect_std_vector_matrix_value() {
  using std::vector;
  using Eigen::MatrixXd;

  MatrixXd d1(3, 2);
  d1 << 1, 2, 3, 4, 5, 6;
  MatrixXd d2(3, 2);
  d2 << 10, 20, 30, 40, 50, 60;
  vector<MatrixXd> d;
  d.push_back(d1);
  d.push_back(d2);
  vector<MatrixXd> fd = F::template apply<vector<MatrixXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].rows(), fd[i].rows());
    EXPECT_EQ(d[i].cols(), fd[i].cols());
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
  }
}


template <typename F>
void expect_std_vector_std_vector_matrix_value() {
  using std::vector;
  using Eigen::MatrixXd;

  MatrixXd d1(3, 2);
  d1 << 1, 2, 3, 4, 5, 6;
  MatrixXd d2(3, 2);
  d2 << 10, 20, 30, 40, 50, 60;
  vector<MatrixXd> d;
  d.push_back(d1);
  d.push_back(d2);

  MatrixXd e1(3, 2);
  e1 << -10, -100, -1000, -10000, -100000, -1000000;
  MatrixXd e2(3, 2);
  e2 << 2, 4, 8, 16, 32, 64;
  vector<MatrixXd> e;
  e.push_back(e1);
  e.push_back(e2);

  vector<vector<MatrixXd> > x;
  x.push_back(d);
  x.push_back(e);

  vector<vector<MatrixXd> > fx = F::template apply<vector<vector<MatrixXd> > >(x);
  EXPECT_EQ(x.size(), fx.size());
  for (size_t i = 0; i < x.size(); ++i) {
    EXPECT_EQ(x[i].size(), fx[i].size());
    for (size_t j = 0; j < x[i].size(); ++j) {
      EXPECT_EQ(x[i][j].rows(), fx[i][j].rows());
      EXPECT_EQ(x[i][j].cols(), fx[i][j].cols());
      for (int k = 0; k < x[i][j].size(); ++k)
        EXPECT_FLOAT_EQ(F::apply_base(x[i][j](k)), fx[i][j](k));
    }
  }
}


// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_values() {
  expect_scalar_value<F>();
  expect_std_vectors_value<F>();
  expect_matrix_value<F>();
  expect_vector_value<F>();
  expect_row_vector_value<F>();
  expect_std_vector_row_vector_value<F>();
  expect_std_vector_vector_value<F>();
  expect_std_vector_matrix_value<F>();
  expect_std_vector_std_vector_matrix_value<F>();
}

#endif

#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

template <typename F>
void expect_scalar_value() {
  using std::vector;  
  vector<double> valid_inputs = F::valid_inputs();
  vector<int> int_valid_inputs = F::int_valid_inputs();
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    double v = F::template apply<double>(input);
    EXPECT_FLOAT_EQ(F::apply_base(double(input)), v);
  }
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    EXPECT_FLOAT_EQ(F::apply_base(valid_inputs[i]), 
                    F::template apply<double>(valid_inputs[i]));
  }
}

template <typename F>
void expect_std_vectors_value() {
  using std::vector;

  vector<double> valid_inputs = F::valid_inputs();

  vector<double> fy = F::template apply<vector<double> >(valid_inputs);
  EXPECT_EQ(valid_inputs.size(), fy.size());
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    EXPECT_FLOAT_EQ(F::apply_base(valid_inputs[i]), fy[i]);
  }

  vector<vector<double> > z;
  z.push_back(valid_inputs);
  z.push_back(valid_inputs);
  vector<vector<double> > fz 
    = F::template apply<vector<vector<double> > >(z);
  EXPECT_EQ(z.size(), fz.size());
  for (size_t i = 0; i < fz.size(); ++i) {
    EXPECT_EQ(z[i].size(), fz[i].size());
    for (size_t j = 0; j < z[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(z[i][j]), fz[i][j]);
  }

  vector<int> int_valid_inputs = F::int_valid_inputs();
  vector<double> u2 = F::template apply<vector<double> >(int_valid_inputs);
  EXPECT_EQ(int_valid_inputs.size(), u2.size());
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    EXPECT_FLOAT_EQ(F::apply_base(int_valid_inputs[i]), u2[i]);
  }
}

template <typename F>
void expect_matrix_value() {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using std::vector;

  vector<double> valid_inputs = F::valid_inputs();
  RowVectorXd eigen_valid_inputs = 
    RowVectorXd::Map(valid_inputs.data(), valid_inputs.size()); 
  size_t num_rows = 3;
  MatrixXd a(num_rows, valid_inputs.size());

  for (size_t i = 0; i < num_rows; i++) {
    a.row(i) << eigen_valid_inputs;
  }

  MatrixXd fa = F::template apply<MatrixXd>(a);
  EXPECT_EQ(a.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(a(i)), fa(i));
  
  MatrixXd fab = F::template apply<MatrixXd>(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(F::apply_base(a(1,1)), fab(0,0));

  vector<MatrixXd> d;
  d.push_back(a);
  d.push_back(a);
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
void expect_vector_value() {
  using Eigen::VectorXd;
  using std::vector;

  std::vector<double> valid_inputs = F::valid_inputs();

  VectorXd b = VectorXd::Map(valid_inputs.data(), valid_inputs.size());
  VectorXd fb = F::template apply<VectorXd>(b);
  EXPECT_EQ(b.size(), fb.size());
  for (int i = 0; i < fb.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(b(i)), fb(i));

  vector<VectorXd> d;
  d.push_back(b);
  d.push_back(b);
  vector<VectorXd> fd = F::template apply<vector<VectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
  }
}

template <typename F>
void expect_row_vector_value() {
  using Eigen::RowVectorXd;
  using std::vector;

  std::vector<double> valid_inputs = F::valid_inputs();

  RowVectorXd c = 
    RowVectorXd::Map(valid_inputs.data(), valid_inputs.size());

  RowVectorXd fc = F::template apply<RowVectorXd>(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i)
    EXPECT_FLOAT_EQ(F::apply_base(c(i)), fc(i));

  vector<RowVectorXd> d;
  d.push_back(c);
  d.push_back(c);

  vector<RowVectorXd> fd = F::template apply<vector<RowVectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j)
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
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
}

#endif

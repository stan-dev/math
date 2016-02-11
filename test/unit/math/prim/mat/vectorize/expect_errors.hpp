#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ERRORS_HPP

#include <stdexcept>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>


template <typename F>
void expect_scalar_error() {
  using std::vector;  

  vector<double> illegal_inputs = F::illegal_inputs();
  vector<int> int_illegal_inputs = F::int_illegal_inputs();

  for (size_t i = 0; i < int_illegal_inputs.size(); ++i) {
    int input = int_illegal_inputs[i];
    EXPECT_THROW(F::template apply<int>(input), std::domain_error);
    EXPECT_THROW(F::apply_base(static_cast<double>(input)), 
      std::domain_error);
  }

  for (size_t i = 0; i < illegal_inputs.size(); ++i) {
    EXPECT_THROW(F::apply_base(illegal_inputs[i]), 
                                               std::domain_error);
  }
}

template <typename F>
void expect_std_vectors_error() {
  using stan::math::foo;
  using std::vector;

  vector<double> illegal_inputs = F::illegal_inputs();

  EXPECT_THROW(F::template apply<vector<double> >(illegal_inputs), 
                                                      std::domain_error);

  vector<vector<double> > z;
  z.push_back(illegal_inputs);
  z.push_back(illegal_inputs);
  EXPECT_THROW(F::template apply<vector<vector<double> > >(z),
                                                     std::domain_error);

  vector<int> int_illegal_inputs = F::int_illegal_inputs();
  EXPECT_THROW(F::template apply<vector<double> >(int_illegal_inputs), 
                                                      std::domain_error);
  vector<vector<int> > z1;
  z1.push_back(int_illegal_inputs);
  z1.push_back(int_illegal_inputs);
  EXPECT_THROW(
    F::template apply<vector<vector<double> > >(z1), std::domain_error);
}

template <typename F>
void expect_matrix_error() {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using std::vector;

  vector<double> illegal_inputs = F::illegal_inputs();
  RowVectorXd eigen_illegal_inputs = 
    RowVectorXd::Map(illegal_inputs.data(), illegal_inputs.size()); 
  size_t num_rows = 3;
  MatrixXd a(num_rows, illegal_inputs.size());

  for (size_t i = 0; i < num_rows; i++) {
    a.row(i) << eigen_illegal_inputs;
  }

  EXPECT_THROW(F::template apply<MatrixXd>(a), std::domain_error);
  
  EXPECT_THROW(F::template apply<MatrixXd>(a.block(1, 1, 1, 1)), 
                                               std::domain_error);

  vector<MatrixXd> d;
  d.push_back(a);
  d.push_back(a);
  EXPECT_THROW(F::template apply<vector<MatrixXd> >(d), std::domain_error);
}

template <typename F>
void expect_vector_error() {
  using Eigen::VectorXd;
  using std::vector;

  std::vector<double> illegal_inputs = F::illegal_inputs();

  VectorXd b = VectorXd::Map(illegal_inputs.data(), illegal_inputs.size());
  EXPECT_THROW(F::template apply<VectorXd>(b), std::domain_error);

  vector<VectorXd> d;
  d.push_back(b);
  d.push_back(b);
  EXPECT_THROW(F::template apply<vector<VectorXd> >(d), std::domain_error);
}

template <typename F>
void expect_row_vector_error() {
  using Eigen::RowVectorXd;
  using std::vector;

  std::vector<double> illegal_inputs = F::illegal_inputs();

  RowVectorXd c = 
    RowVectorXd::Map(illegal_inputs.data(), illegal_inputs.size());

  EXPECT_THROW(F::template apply<RowVectorXd>(c), std::domain_error);

  vector<RowVectorXd> d;
  d.push_back(c);
  d.push_back(c);

  EXPECT_THROW(F::template apply<vector<RowVectorXd> >(d), 
                                                  std::domain_error);
}

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_errors() {
  expect_scalar_error<F>();
  expect_std_vectors_error<F>();
  expect_matrix_error<F>();
  expect_vector_error<F>();
  expect_row_vector_error<F>();
}

#endif

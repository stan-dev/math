#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ERRORS_HPP

#include <stdexcept>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <test/unit/math/fwd/mat/vectorize/expect_match_return_t.hpp>
#include <stan/math/fwd/core.hpp>

template <typename F>
void expect_std_vectors_error() {
  using std::vector;
  using stan::math::fvar;

  vector<double> illegal_inputs = F::illegal_inputs();
  vector<fvar<double> > y;

  for (size_t i = 0; i < illegal_inputs.size(); ++i) {
    EXPECT_THROW(F::template apply<fvar<double> >(illegal_inputs[i]), 
                   std::domain_error);
    y.push_back(illegal_inputs[i]);
  }

  EXPECT_THROW(F::template apply<vector<fvar<double> > >(y), std::domain_error);

  vector<vector<fvar<double> > > z;
  z.push_back(y);
  z.push_back(y);
  EXPECT_THROW(F::template apply<vector<vector<fvar<double> > > >(z),
                                                     std::domain_error);
}

template <typename F>
void expect_matrix_error() {
  using std::vector;
  using stan::math::fvar;
  typedef 
    Eigen::Matrix<fvar<double> , Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;

  vector<double> illegal_inputs = F::illegal_inputs();
  size_t num_rows = 3;
  MatrixXvar a(num_rows, illegal_inputs.size());

  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < illegal_inputs.size(); ++j)
      a(i, j) = illegal_inputs[j];
  }

  EXPECT_THROW(F::template apply<MatrixXvar>(a), std::domain_error);
  
  EXPECT_THROW(F::template apply<MatrixXvar>(a.block(1, 1, 1, 1)), 
                                               std::domain_error);

  vector<MatrixXvar> d;
  d.push_back(a);
  d.push_back(a);
  EXPECT_THROW(F::template apply<vector<MatrixXvar> >(d), 
                                               std::domain_error);
}

template <typename F>
void expect_vector_error() {
  using std::vector;
  using stan::math::fvar;
  typedef Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> VectorXvar;

  std::vector<double> illegal_inputs = F::illegal_inputs();

  VectorXvar b = VectorXvar(illegal_inputs.size());
  for (size_t i = 0; i < illegal_inputs.size(); ++i) 
    b(i) = illegal_inputs[i];
  EXPECT_THROW(F::template apply<VectorXvar>(b), std::domain_error);

  vector<VectorXvar> d;
  d.push_back(b);
  d.push_back(b);
  EXPECT_THROW(F::template apply<vector<VectorXvar> >(d), 
                                                 std::domain_error);
}

template <typename F>
void expect_row_vector_error() {
  using std::vector;
  using stan::math::fvar;
  typedef Eigen::Matrix<fvar<double> , 1, Eigen::Dynamic> RowVectorXvar;

  std::vector<double> illegal_inputs = F::illegal_inputs();

  RowVectorXvar c = RowVectorXvar(illegal_inputs.size());
  for (size_t i = 0; i < illegal_inputs.size(); ++i) 
    c(i) = illegal_inputs[i];
  EXPECT_THROW(F::template apply<RowVectorXvar>(c), std::domain_error);

  vector<RowVectorXvar> d;
  d.push_back(c);
  d.push_back(c);

  EXPECT_THROW(F::template apply<vector<RowVectorXvar> >(d), 
                                                  std::domain_error);
}

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_errors() {
  expect_std_vectors_error<F>();
  expect_matrix_error<F>();
  expect_vector_error<F>();
  expect_row_vector_error<F>();
}

#endif

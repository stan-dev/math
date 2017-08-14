#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP

#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

template <typename F, typename matrix_t, typename matrix_d>
void expect_binary_matrix_size_err_throw(matrix_t template_tm,
                                           matrix_d template_dm) {
  using std::vector;
  using Eigen::MatrixXd;

  matrix_t badsize_tm1 = build_template_matrix(template_tm, 4, 5);
  matrix_t badsize_tm2 = build_template_matrix(template_tm, 7, 5);
  matrix_t badsize_tm3 = build_template_matrix(template_tm, 4, 9);
  matrix_t badsize_tm4 = build_template_matrix(template_tm, 3, 6);
  matrix_t badsize_tm5 = build_template_matrix(template_tm, 5, 4);
  matrix_d badsize_dm1 = build_template_matrix(template_dm, 6, 13);
  matrix_d badsize_dm2 = build_template_matrix(template_dm, 4, 8);
  matrix_d badsize_dm3 = build_template_matrix(template_dm, 12, 5);
  matrix_d badsize_dm4 = build_template_matrix(template_dm, 5, 4);

  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm1, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm2),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm2, badsize_tm1),
  std::invalid_argument);

  if (badsize_tm1.rows() > 1 && badsize_tm1.cols() > 1) {
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm2),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm2, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm3),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm3, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm4),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm4, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm3),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm3, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm4),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm4, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm5),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm5, badsize_tm1),
                 std::invalid_argument);
  }
}

template <typename F, typename matrix_t, typename matrix_d>
void expect_binary_matrix_value_err_throw(matrix_t template_tm,
                                          matrix_d template_dm) {
  using std::vector;

  vector<double> invalid_inputs1 = F::invalid_inputs1();
  if (invalid_inputs1.size() == 0) return;
  vector<double> invalid_inputs2 = F::invalid_inputs2();
  vector<int> int_invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> int_invalid_inputs2 = F::int_invalid_inputs2();
  matrix_d a1 = build_template_matrix(template_dm, 
                                      invalid_inputs1.size(), 3);
  matrix_t b1 = build_template_matrix(template_tm, 
                                      invalid_inputs1.size(), 3);
  for (int i = 0; i < a1.size(); ++i) {
    a1(i) = invalid_inputs1[(i % invalid_inputs1.size())];
    b1(i) = invalid_inputs1[(i % invalid_inputs1.size())];
  }
  matrix_d a2 = build_template_matrix(template_dm, 
                                      invalid_inputs2.size(), 3);
  matrix_t b2 = build_template_matrix(template_tm, 
                                      invalid_inputs2.size(), 3);
  for (int i = 0; i < a2.size(); ++i) {
    a2(i) = invalid_inputs2[(i % invalid_inputs2.size())];
    b2(i) = invalid_inputs2[(i % invalid_inputs2.size())];
  }

  EXPECT_THROW(F::template apply<matrix_t>(a1, b2), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1, a2), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1, b2), std::domain_error);
  if (a1.rows() > 1 && a2.cols() > 1) {
    EXPECT_THROW(F::template apply<matrix_t>(b1.block(1, 1, 1, 1),
                                             a2.block(1, 1, 1, 1)),
                 std::domain_error);
    EXPECT_THROW(F::template apply<matrix_t>(a1.block(1, 1, 1, 1),
                                             b2.block(1, 1, 1, 1)),
                 std::domain_error);
    EXPECT_THROW(F::template apply<matrix_t>(b1.block(1, 1, 1, 1),
                                             b2.block(1, 1, 1, 1)),
                 std::domain_error);
  }

  vector<matrix_d> d1;
  d1.push_back(a1);
  d1.push_back(a1);
  vector<matrix_d> d2;
  d2.push_back(a2);
  d2.push_back(a2);
  vector<matrix_t> e1;
  e1.push_back(b1);
  e1.push_back(b1);
  vector<matrix_t> e2;
  e2.push_back(b2);
  e2.push_back(b2);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(e1, d2),
               std::domain_error);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(d1, e2),
               std::domain_error);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(e1, e2),
               std::domain_error);
}


template <typename F, typename T, int R, int C>
void expect_binary_matrix_error(Eigen::Matrix<double, R, C> template_m) {
  Eigen::Matrix<T, R, C> template_tm;
  expect_binary_matrix_size_err_throw<F>(template_tm, template_m);
  expect_binary_matrix_value_err_throw<F>(template_tm, template_m);
}

#endif

#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP

#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_matrix_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_std_vector_matrix_err_throw.hpp>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

template <typename F, typename V>
void expect_row_vector_error() {
  using std::vector;
  using Eigen::RowRowVectorXd;
  typedef Eigen::Matrix<V, 1, Eigen::Dynamic> row_row_vector_t;
  row_vector_t badsize_tm1(4);
  row_vector_t badsize_tm2(9);
  RowVectorXd badsize_dm(11);

  EXPECT_THROW(F::template apply<row_vector_t>(badsize_tm1, badsize_dm),
               std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(badsize_dm, badsize_tm1),
               std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(badsize_tm1, badsize_tm2),
               std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(badsize_tm2, badsize_tm1),
               std::domain_error);

  vector<double> invalid_inputs = F::invalid_inputs();
  if (invalid_inputs.size() == 0) return;
  vector<int> int_invalid_inputs = F::int_invalid_inputs();
  vector<T> y(invalid_inputs.begin(), invalid_inputs.end());
  RowVectorXd a(invalid_inputs.size());
  row_vector_t b(invalid_inputs.size());
  for (int i = 0; i < a.size(); ++i) {
      a(i) = invalid_inputs[i];
      b(i) = invalid_inputs[i];
    }
  }
  expect_binary_scalar_matrix_err_throw<F>(y, a);
  expect_binary_scalar_matrix_err_throw<F>(y, b);
  expect_binary_scalar_matrix_err_throw<F>(int_invalid_inputs, b);
  expect_binary_scalar_matrix_err_throw<F>(invalid_inputs, b);

  EXPECT_THROW(F::template apply<row_vector_t>(a, b), std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(b, a), std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(b, b), std::domain_error);
  EXPECT_THROW(F::template apply<row_vector_t>(b.block(1, 1, 1, 1),
                                               b.block(1, 1, 1, 1)),
               std::domain_error);

  vector<RowRowVectorXd> d;
  d.push_back(a);
  d.push_back(a);
  vector<row_vector_t> e;
  e.push_back(b);
  e.push_back(b);
  expect_binary_scalar_std_vector_matrix_err_throw<F>(y, d);
  expect_binary_scalar_std_vector_matrix_err_throw<F>(y, e);
  expect_binary_scalar_std_vector_matrix_err_throw<F>(
    int_invalid_inputs, e);
  expect_binary_scalar_std_vector_matrix_err_throw<F>(invalid_inputs, e);
  EXPECT_THROW(F::template apply<row_vector_t>(d, e), std::domain_error);
}
#endif

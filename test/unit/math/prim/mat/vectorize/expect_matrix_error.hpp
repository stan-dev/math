#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_MATRIX_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_MATRIX_ERROR_HPP

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <exception>
#include <vector>

template <typename F, typename T>
void expect_matrix_error() {
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
  vector<double> invalid_inputs = F::invalid_inputs();
  if (invalid_inputs.size() == 0) return;
  matrix_t a(3, invalid_inputs.size());
  for (int i = 0; i < a.rows(); ++i)
    for (int j = 0; j < a.cols(); ++j)
      a(i, j) = invalid_inputs[j];
  EXPECT_THROW(F::template apply<matrix_t>(a), std::exception);
  EXPECT_THROW(F::template apply<matrix_t>(a.block(1, 1, 1, 1)), 
               std::exception);

  vector<matrix_t> d;
  d.push_back(a);
  d.push_back(a);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(d), 
               std::exception);
}

#endif

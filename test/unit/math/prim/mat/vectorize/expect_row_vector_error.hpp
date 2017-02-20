#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP

#include <Eigen/Dense>
#include <exception>
#include <gtest/gtest.h>

template <typename F, typename V>
void expect_row_vector_error() {
  using std::vector;
  typedef Eigen::Matrix<V, 1, Eigen::Dynamic> row_vector_t;
  std::vector<double> invalid_inputs = F::invalid_inputs();
  if (invalid_inputs.size() == 0) return;
  row_vector_t c = row_vector_t(invalid_inputs.size());
  for (size_t i = 0; i < invalid_inputs.size(); ++i) 
    c(i) = invalid_inputs[i];
  EXPECT_THROW(F::template apply<row_vector_t>(c), std::exception);

  vector<row_vector_t> d;
  d.push_back(c);
  d.push_back(c);

  EXPECT_THROW(F::template apply<vector<row_vector_t> >(d), 
               std::exception);
}
#endif

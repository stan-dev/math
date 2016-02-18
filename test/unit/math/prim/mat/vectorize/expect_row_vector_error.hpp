#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_ROW_VECTOR_ERROR_HPP

#include <Eigen/Dense>
#include <stdexcept>
#include <gtest/gtest.h>

template <typename F, typename V>
void expect_row_vector_error() {
  using std::vector;
  typedef Eigen::Matrix<V, 1, Eigen::Dynamic> row_vector_t;

  std::vector<double> illegal_inputs = F::illegal_inputs();

  row_vector_t c = row_vector_t(illegal_inputs.size());
  for (size_t i = 0; i < illegal_inputs.size(); ++i) 
    c(i) = illegal_inputs[i];
  EXPECT_THROW(F::template apply<row_vector_t>(c), std::domain_error);

  vector<row_vector_t> d;
  d.push_back(c);
  d.push_back(c);

  EXPECT_THROW(F::template apply<vector<row_vector_t> >(d), 
               std::domain_error);
}
#endif

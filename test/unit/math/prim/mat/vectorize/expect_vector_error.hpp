#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VECTOR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VECTOR_HPP

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <exception>
#include <vector>

template <typename F, typename T>
void expect_vector_error() {
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_t;
  vector<double> invalid_inputs = F::invalid_inputs();
  if (invalid_inputs.size() == 0) return;
  vector_t b = vector_t(invalid_inputs.size());
  for (size_t i = 0; i < invalid_inputs.size(); ++i) 
    b(i) = invalid_inputs[i];
  EXPECT_THROW(F::template apply<vector_t>(b), std::exception);

  vector<vector_t> d;
  d.push_back(b);
  d.push_back(b);
  EXPECT_THROW(F::template apply<vector<vector_t> >(d), 
               std::exception);
}

#endif

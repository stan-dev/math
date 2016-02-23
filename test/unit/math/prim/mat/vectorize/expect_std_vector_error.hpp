#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_STD_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_STD_VECTOR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F, typename T>
void expect_std_vector_error() {
  using std::vector;
  vector<double> invalid_inputs = F::invalid_inputs();
  vector<T> y(invalid_inputs.begin(), invalid_inputs.end());
  EXPECT_THROW(F::template apply<vector<T> >(y), std::domain_error);

  vector<vector<T> > z;
  z.push_back(y);
  z.push_back(y);
  EXPECT_THROW(F::template apply<vector<vector<T> > >(z),
               std::domain_error);
}

#endif

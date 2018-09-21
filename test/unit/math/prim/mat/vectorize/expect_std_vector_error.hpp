#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_STD_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_STD_VECTOR_ERROR_HPP

#include <gtest/gtest.h>
#include <exception>
#include <vector>

template <typename F, typename T>
void expect_std_vector_error() {
  using std::vector;
  vector<double> invalid_inputs = F::invalid_inputs();
  if (invalid_inputs.size() == 0)
    return;
  vector<T> y(invalid_inputs.begin(), invalid_inputs.end());
  EXPECT_THROW(F::template apply<vector<T> >(y), std::exception);

  vector<vector<T> > z;
  z.push_back(y);
  z.push_back(y);
  EXPECT_THROW(F::template apply<vector<vector<T> > >(z), std::exception);
}

#endif

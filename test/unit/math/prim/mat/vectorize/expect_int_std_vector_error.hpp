#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_STD_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_STD_VECTOR_ERROR_HPP

#include <gtest/gtest.h>
#include <exception>
#include <vector>

template <typename F>
void expect_int_std_vector_error() {
  using std::vector;
  vector<int> invalid_inputs = F::int_invalid_inputs();
  if (invalid_inputs.size() == 0) return;
  EXPECT_THROW(F::template apply<vector<double> >(invalid_inputs), 
               std::exception);

  vector<vector<int> > z;
  z.push_back(invalid_inputs);
  z.push_back(invalid_inputs);
  EXPECT_THROW(F::template apply<vector<vector<double> > >(z),
               std::exception);
}

#endif

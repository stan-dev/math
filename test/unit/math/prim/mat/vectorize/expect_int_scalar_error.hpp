#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_SCALAR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_SCALAR_ERROR_HPP

#include <gtest/gtest.h>
#include <exception>
#include <vector>

template <typename F>
void expect_int_scalar_error() {
  using std::vector;
  vector<int> int_invalid_inputs = F::int_invalid_inputs();
  for (size_t i = 0; i < int_invalid_inputs.size(); ++i) {
    int input = int_invalid_inputs[i];
    EXPECT_THROW(F::template apply<double>(input), std::exception);
  }
}

#endif

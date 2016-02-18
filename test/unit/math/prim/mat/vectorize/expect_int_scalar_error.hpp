#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_SCALAR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_INT_SCALAR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F>
void expect_int_scalar_error() {
  using std::vector;  
  vector<int> int_illegal_inputs = F::int_illegal_inputs();
  for (size_t i = 0; i < int_illegal_inputs.size(); ++i) {
    int input = int_illegal_inputs[i];
    EXPECT_THROW(F::template apply<double>(input), std::domain_error);
  }
}

#endif

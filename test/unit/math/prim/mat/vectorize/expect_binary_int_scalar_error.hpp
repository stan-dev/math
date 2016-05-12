#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_INT_SCALAR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_INT_SCALAR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F>
void expect_binary_int_scalar_error() {
  using std::vector;  
  vector<int> int_invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> int_invalid_inputs2 = F::int_invalid_inputs2();
  for (size_t i = 0; i < int_invalid_inputs1.size(); ++i) {
    EXPECT_THROW(F::template apply<double>(int_invalid_inputs1[i], 
    int_invalid_inputs2[i]), std::domain_error);
  }
}

#endif

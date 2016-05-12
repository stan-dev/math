#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F, typename T>
void expect_binary_scalar_error() {
  using std::vector;  
  vector<int> int_invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> int_invalid_inputs2 = F::int_invalid_inputs2();
  vector<double> invalid_inputs1 = F::invalid_inputs1();
  vector<double> invalid_inputs2 = F::invalid_inputs2();
  for (size_t i = 0; i < invalid_inputs1.size(); ++i) {
    T input1 = invalid_inputs1[i];
    T input2 = invalid_inputs2[i];
    EXPECT_THROW(F::template apply<T>(input1, int_invalid_inputs2[i]), 
    std::domain_error);
    EXPECT_THROW(F::template apply<T>(int_invalid_inputs1[i], input2), 
    std::domain_error);
    EXPECT_THROW(F::template apply<T>(input1, invalid_inputs2[i]), 
    std::domain_error);
    EXPECT_THROW(F::template apply<T>(invalid_inputs1[i], input2), 
    std::domain_error);
    EXPECT_THROW(F::template apply<T>(input1, input2), std::domain_error);
  }
}

#endif

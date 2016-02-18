#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_SCALAR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_SCALAR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F, typename T>
void expect_scalar_error() {
  using std::vector;  
  vector<double> illegal_inputs = F::illegal_inputs();
  for (size_t i = 0; i < illegal_inputs.size(); ++i) {
    T input = illegal_inputs[i];
    EXPECT_THROW(F::template apply<T>(input), std::domain_error);
  }
}

#endif
